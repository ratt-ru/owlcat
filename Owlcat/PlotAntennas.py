# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os

import casacore.measures
import casacore.quanta as qa
import casacore.measures
import numpy as np

from casacore.tables import table
from collections import namedtuple

from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, PreText
from bokeh.io import output_file, save


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)-15s %(filename)s %(levelname)s %(message)s")

logging.disable(logging.CRITICAL)


def wgs84_to_ecef(lon, lat, alt):
    """ 
    Convert wgs84(latitude (deg), longitude(deg), elevation(deg)) to
    Earth Centred Earth fixed coordinates (X(m), Y(m), Z(m)).

    Coordinates in the ITRF format should be converted to WGS84 coordinate system for consistency.
    This function is an implementation based on the following reference:

    https://docs.hisparc.nl/coordinates/HiSPARC_coordinates.pdf

    Parameters
    ----------:
    lat: :obj:`float`
        Latitude in degrees
    lon: :obj:`float`
        Longitude in degrees
    alt: :obj:`float`
        Altitude in metres

    Returns
    -------
    X, Y, Z:  ECEF coordinates in metres
    """

    logging.debug("Converting WGS84 to ECF")

    # set up earth's shape ellipsoid approximation
    # semi major axis
    a = 6378137.0

    # flattening
    f = 1 / 298.257223563

    # semi-minor axis
    b = a - a * f

    # eccentricity
    e = np.sqrt((2 * f) - (f**2))

    # Normal: Distance between a location on the ellipsoid and the
    # intersection of its nromal and the ellipsoid's z-axis
    N = a / np.sqrt(1 - e**2 * np.sin(lat)**2)

    # altitude
    h = alt

    X = (N + h) * np.cos(lat) * np.cos(lon)
    Y = (N + h) * np.cos(lat) * np.sin(lon)
    Z = (((b**2 / a**2) * N) + h) * np.sin(lat)

    return X, Y, Z


def ecef_to_enu(x, y, z, r_lon, r_lat, r_alt):
    """
    Convert ECEF coordinates to East, North, Up (m). Because the ENU is a 
    local coordinate system (to the observer/ antenna in this case)
    an reference position / origin is required from which calculate the 
    relative positions.

    Parameters
    ----------
    x: :obj:`float`
        X in metres
    y: :obj:`float`
        Y in metres
    z: :obj:`float`
        Z in metres
    rx: :obj:`float`
        Reference latitude in degrees
    ry: :obj:`float`
        Reference longitude in degrees

    Returns
    -------
    ee, en, eu: :obj:`float`
        East, North Up positions in metres
    """

    logging.debug("Converting ECEF to ENU")

    # convert reference coordinates to ECEF
    rx, ry, rz = wgs84_to_ecef(r_lon, r_lat, r_alt)

    dx = x - rx
    dy = y - ry
    dz = z - rz

    ee = (-np.sin(r_lon) * dx) \
        + (np.cos(r_lon) * dy) \
        + (0 * dz)

    en = (-np.sin(r_lat) * np.cos(r_lon) * dx) \
        + (-np.sin(r_lat) * np.sin(r_lon) * dy) \
        + (np.cos(r_lat) * dz)

    eu = (np.cos(r_lat) * np.cos(r_lon) * dx) \
        + (np.cos(r_lat) * np.sin(r_lon) * dy) \
        + (np.sin(r_lat) * dz)

    return ee, en, eu


def wgs84_to_enu(lon, lat, alt, r_lon, r_lat, r_alt):
    """Convert WSG84 / ITRF to ENU"""
    X, Y, Z = wgs84_to_ecef(lon, lat, alt)
    e, n, u = ecef_to_enu(X, Y, Z, r_lon, r_lat, r_alt)

    return e, n, u


def get_antenna_data(ms_name):
    """Return antenna ITRF positions and other data"""

    logging.debug("Getting antenna data")

    ms = table(ms_name, ack=False)

    # get antenna related information
    ant_sub = table(ms.getkeyword("ANTENNA"), ack=False)
    names = ant_sub.getcol("NAME")
    stations = ant_sub.getcol("STATION")
    offsets = ant_sub.getcol("OFFSET")
    positions = ant_sub.getcol("POSITION")
    indices = ant_sub.rownumbers()
    ant_sub.close()

    # get the name of the observatory / telescope
    observ_sub = table(ms.getkeyword("OBSERVATION"), ack=False)
    obs = observ_sub.getcell("TELESCOPE_NAME", 0)
    observ_sub.close()

    ms.close()

    logging.debug(f"Telescope name: {obs}")
    logging.debug(f"Found: {len(names)} antennas")

    meta = namedtuple("antenna",
                      "names, stations, positions, offsets, telescope, indices")
    antenna = meta(names=names, stations=stations, positions=positions,
                   offsets=offsets, telescope=obs, indices=indices)

    return antenna


def get_antenna_coords(positions, offsets, obs):
    """Return antenna COFA coordinates, and antenna coordinates

    Parameters
    ----------
    positions: :obj:`np.array`
        Antenna ITRF positions
    offset: :obj:`np.array`
        Antenna offset positions
    obs: :obj:`str`
        Telescope name or observatory name
    """

    logging.debug("Getting antenna coordinates")

    me = casacore.measures.measures()

    x, y, z = positions.T.tolist()
    x_off, y_off, z_off = offsets.T.tolist()

    # make ITRF positions casacore quantities so that we can process using
    # casacore measurements
    xq = qa.quantity(x, 'm')
    yq = qa.quantity(y, 'm')
    zq = qa.quantity(z, 'm')

    xq_off = qa.quantity(x_off, 'm')
    yq_off = qa.quantity(y_off, 'm')
    zq_off = qa.quantity(z_off, 'm')

    # make antenna positions into casacore measurements
    baseline_me = me.position('itrf', v0=xq, v1=yq, v2=zq,
                              off=[xq_off, yq_off, zq_off])

    # get the list of observatories
    obs_list = me.get_observatories()

    if obs in obs_list:
        # Get a (position) measure Centre of Array (cofa) for the given
        # osbervatory. This is in WSG84 reference frame.

        obs_cofa = me.observatory(obs)

        # get the cofa's reference frame
        obs_cofa = me.measure(obs_cofa, "WGS84")

        # make this position the reference frame of the entire measurement
        me.doframe(obs_cofa)
    else:
        logging.error(f"COFA of Observatory: {obs} not found")

    # latitude, longitude and elevation in degree, minute second string form
    # get this before changing reference frame
    str_lon_lat_el = me.get_value(baseline_me)
    str_lon_lat_el = [_.formatted("dd.mm.ss.t..").strip(
        "[]").split(",") for _ in str_lon_lat_el]

    # change the baseline measurement to WGS84
    baseline_me = me.measure(baseline_me, "WGS84")

    # get the wsg84 coordinates in radians for the COFA
    cofa_coords = me.get_value(obs_cofa)
    cofa_coords = [_.get_value() for _ in cofa_coords]

    # for the antennas
    bl_coords = me.get_value(baseline_me)
    bl_coords = [_.get_value() for _ in bl_coords]

    rf_coords = namedtuple("coords", "cofa, ant_rad, ant_str")
    coords = rf_coords(cofa=cofa_coords, ant_rad=bl_coords,
                       ant_str=str_lon_lat_el)

    return coords


def plot_antennas(source, ms_name, nants):
    """Make bokeh plot of the antennas

    Parameters
    ----------
    source: :obj:`ColumnDataSource`
        Data source containing the plot's data

    Return
    ------
    out_layout: Bokeh layout object containing the plots
    """

    tooltips = [
        ("Antenna", "@name"),
        ("Station", "@station"),
        ("Antenna index", "@indices"),
        ("Lon, Lat [dms]", "@lon, @lat"),
        ("E,N,U [m]", "(@e{0.000}, @n{0.000}, @u{0.000})"),
        ("ITRF (x, y, z) [m]", "(@x{0.000}, @y{0.000}, @z{0.000})"),
    ]

    logging.debug("Initialising figure")

    r_max = 1.1 * max(max(source["e"]), max(source["n"]))
    r_min = 1.1 * min(min(source["e"]), min(source["n"]))
    # ranges = (r_min, r_max)
    ranges = None

    fig = figure(plot_width=720, plot_height=720, sizing_mode="scale_height",
                 x_axis_label="East [m]", y_axis_label="North [m]",
                 title="Offset From Array Centre", toolbar_location="above",
                 tooltips=tooltips, x_range=ranges, y_range=ranges)

    fig.update(outline_line_width=2, outline_line_color="#017afe",
               outline_line_alpha=0.4)

    fig.title.update(align="center", text_font_size="15pt")

    fig.xaxis.update(axis_label_text_font="monospace",
                     axis_label_text_font_size="10pt",
                     axis_label_text_font_style="normal")
    fig.yaxis.update(axis_label_text_font="monospace",
                     axis_label_text_font_size="10pt",
                     axis_label_text_font_style="normal")

    color = source.pop("color")

    source = ColumnDataSource(data=source)
    p = fig.scatter("e", "n", line_width=2, line_color=color,
                    marker="circle", fill_color=None, size=30, fill_alpha=0.8,
                    source=source)
    fig.text(x="e", y="n", text="name", text_align="center",
             text_baseline="middle", text_font_size="7pt", source=source)

    fig.hover.renderers = [p]

    infos = f"MS      : {ms_name}\n"
    infos += f"Antennas: {nants}"
    pre = PreText(text=infos, sizing_mode="stretch_width")

    out_layout = column(children=[pre, fig], sizing_mode="stretch_both")

    return out_layout


def main(args):

    ms_name = args.ms_name

    if args.o_file:
        o_file = args.o_file
    else:
        o_file = os.path.basename(os.path.abspath(ms_name))
        o_file += "_antennas.html"

    output_file(o_file)

    # get colors
    if args.cmap:
        cmap = args.cmap
    else:
        cmap = "#FC8103"

    # get information from antenna data table
    ant = get_antenna_data(ms_name)

    # colors = all_palettes[cmap][max(all_palettes[cmap].keys())]
    # colors = linear_palette(colors, len(ant.names))

    # transform antenna coordinates
    coords = get_antenna_coords(ant.positions, ant.offsets, ant.telescope)

    cofa_lon, cofa_lat, cofa_alt = coords.cofa

    # antenna coords in lon / lat in radians
    ant_lon, ant_lat, ant_alt = coords.ant_rad

    # antenna coords in lon/ lat in string format
    str_lon, str_lat, str_el = coords.ant_str

    # ITRF positions
    x, y, z = ant.positions.T.tolist()

    e, n, u = wgs84_to_enu(ant_lon, ant_lat, ant_alt, cofa_lon, cofa_lat,
                           cofa_alt)

    source = dict(e=e, n=n, u=u,
                  lon=str_lon, lat=str_lat, el=str_el,
                  x=x, y=y, z=z,
                  name=ant.names, station=ant.stations, indices=ant.indices,
                  color=cmap)

    plot = plot_antennas(source, ms_name, len(ant.names))

    save(plot)

    print(f"Antenna plot at: {o_file}")
