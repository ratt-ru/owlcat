# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import sys
import logging
import os

import casacore.measures
import casacore.quanta as qa
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
        Latitude in radians
    lon: :obj:`float`
        Longitude in radians
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
    r_lon: :obj:`float`
        Reference longitude in radians
    r_lat: :obj:`float`
        Reference latitude in radians
    r_alt: :obj:`float`
        Reference altitude in metres

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

    with table(ms_name, ack=False) as ms: 
        if "OBSERVATION" in ms.keywordnames():
            with table("::".join((ms_name, "OBSERVATION")), ack=False) as observ_sub:
                obs = observ_sub.getcell("TELESCOPE_NAME", 0)
        else:
            obs = "unknown"

        if "ANTENNA" in ms.keywordnames():
            ant_tab = ms.getkeyword("ANTENNA").split()[-1]
        else:
            ant_tab = ms_name

        with table(ant_tab, ack=False) as ant_sub:
            names = ant_sub.getcol("NAME")
            stations = ant_sub.getcol("STATION")
            offsets = ant_sub.getcol("OFFSET")
            positions = ant_sub.getcol("POSITION")
            indices = ant_sub.rownumbers()

    logging.debug(f"Telescope name: {obs}")
    logging.debug(f"Found: {len(names)} antennas")

    meta = namedtuple("antenna",
                      "names, stations, positions, offsets, telescope, indices")
    antenna = meta(names=names, stations=stations, positions=positions,
                   offsets=offsets, telescope=obs, indices=indices)

    return antenna


def get_antenna_coords(positions, obs, offsets=None, cofa=None):
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

    quantify = lambda coord, unit: qa.quantity(coord, unit)

    # make ITRF positions casacore quantities so that we can process using
    # casacore measurements

    if casacore.measures.is_measure(positions):
        baseline_me = positions

    else:
        xq, yq, zq = list(map(quantify, positions.T, ["m", "m", "m"]))

        if offsets is not None:
            offsets = list(map(quantify, offsets.T, ["m", "m", "m"]))

        # make antenna positions into casacore measurements
        baseline_me = me.position('itrf', v0=xq, v1=yq, v2=zq,
                                  off=offsets)

    if cofa:
        if casacore.measures.is_measure(cofa):
            obs_cofa = cofa
        else:
            cofa = list(map(quantify, [float(_) for _ in cofa.split(",")],
                            ["deg", "deg", "m"]))
            obs_cofa = me.position("WGS84", *cofa)
    else:
        # get the list of observatories
        obs_list = me.get_observatories()

        if obs in obs_list:
            # Get a (position) measure Centre of Array (cofa) for the given
            # osbervatory. This is in WSG84 reference frame.
            obs_cofa = me.observatory(obs)
        else:
            # first antenna position
            i_ant = positions[0]
            ixq, iyq, izq = list(map(quantify, i_ant, ["m", "m", "m"]))
            obs_cofa = me.position('itrf', v0=ixq, v1=iyq, v2=izq,
                                   off=offsets)
            logging.debug(f"Centre of Array of Observatory: {obs} not found.")
            logging.debug("Defaulting to first antenna of the array")

    # get the cofa's reference frame
    obs_cofa = me.measure(obs_cofa, "WGS84")

    logging.debug("COFA: " +
                  ", ".join([_.formatted("dd.mm.ss.t..")
                             for _ in me.get_value(obs_cofa)]))

    # make this position the reference frame of the entire measurement
    me.doframe(obs_cofa)

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


def read_underlay_file(in_file):
    """Read file containing complete antenna layout for plotting.

    An underlay file is a JSON file containing the coordinates of all of a
    telescope's antennas. This includes the antennas that were not used for
    the observation. By default, this script will only plot the antennas that
    were used in the observation i.e. in the current MS. If the underlay file
    is provided, antennas used in the observation will be shown in colour, in
    -30.54417addition to the rest that were inactive, which will be denoted
    by the colour grey.

    Parameters
    ----------
    in_file: :obj:`str`
        Input JSON file containing the telescope's antenna coordinates
    """
    # allowed units
    allowed_us = {"deg", "m", "rad"}

    # allowed rfs
    allowed_rfs = {"WGS84", "ITRF"}

    with open(in_file, "r") as rf:
        data = json.load(rf)

    for spec in ["cofa", "antenna"]:
        rf = [data[spec]["rf"]]
        unit = data[spec]["units"]

        if not set(rf) <= allowed_rfs:
            print(f"Invalid rf {rf} for {spec}")
            sys.exit(-1)

        if not set(unit) <= allowed_us:
            print(f"Invalid unit {unit} for {spec}")
            sys.exit(-1)

    return data


def get_underlay_data(in_file):
    """Parse underlay file data and create position measures for the
    Centre of Array (COFA) and the antenna coordinates
    """
    data = read_underlay_file(in_file)

    obs = data["obs"]

    cofa_units = data["cofa"]["units"]
    cofa_rf = data["cofa"]["rf"]
    cofa_coords = data["cofa"]["coords"]

    ant_units = data["antenna"]["units"]
    ant_rf = data["antenna"]["rf"]
    ant_pos = data["antenna"]["coords"]
    names = list(ant_pos.keys())
    ant_coords = np.array(list(ant_pos.values()))

    quantify = lambda coord, unit: qa.quantity(coord, unit)
    me = casacore.measures.measures()

    c_lon, c_lat, c_el = list(map(quantify, cofa_coords, cofa_units))
    cofa_me = me.position(rf=cofa_rf, v0=c_lon, v1=c_lat, v2=c_el)

    a_lon, a_lat, a_el = list(map(quantify, ant_coords.T, ant_units))
    ant_me = me.position(rf=ant_rf, v0=a_lon, v1=a_lat, v2=a_el)

    underlay = namedtuple("underlay",
                          "cofa, names, positions, offsets, telescope")
    u_data = underlay(names=names, positions=ant_me, cofa=cofa_me,
                      offsets=None, telescope=obs)

    return u_data


def plot_antennas(source, ms_name):
    """Make bokeh plot of the antennas

    Parameters
    ----------
    source: :obj:`ColumnDataSource`
        Data source containing the plot's data
    ms_name: :obj:`str`
        Name of the current MS

    Return
    ------
    out_layout: Bokeh layout object containing the plots
    """

    tooltips = [
        ("Antenna", "@name"),
        ("Station", "@station"),
        ("Index in MS", "@indices"),
        ("Index in full array", "@rindices"),
        ("Lon, Lat [dms]", "@lon, @lat"),
        ("E,N,U [m]", "(@e{0.000}, @n{0.000}, @u{0.000})"),
        ("ITRF (x, y, z) [m]", "(@x{0.000}, @y{0.000}, @z{0.000})"),
    ]

    logging.debug("Initialising figure")

    nants = len(source["name"])
    obs = source.pop("obs")
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

    # create reference data source object
    r_source = {}
    for refs in ["re", "rn", "ru", "rname"]:
        r_source[refs] = source.pop(refs)

    r_source = ColumnDataSource(data=r_source)

    source = ColumnDataSource(data=source)

    fig.circle(x="re", y="rn", source=r_source, color="grey", fill_color=None,
               size=30)
    p = fig.scatter("e", "n", line_width=2, line_color=color, marker="circle",
                    fill_color=None, size=30, source=source)
    fig.text(x="re", y="rn", text="rname", text_align="center",
             text_baseline="middle", text_font_size="7pt", source=r_source)

    fig.hover.renderers = [p]

    infos = f"MS       : {ms_name}\n"
    infos += f"Telescope: {obs}\n"
    infos += f"Antennas : {nants}"
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

    # get colors or set default active colour
    if args.cmap:
        cmap = args.cmap
    else:
        cmap = "#FC8103"

    # get information from antenna data table
    ant = get_antenna_data(ms_name)

    # gather antenna coordinates
    coords = get_antenna_coords(ant.positions, ant.telescope,
                                offsets=ant.offsets, cofa=args.cofa)

    # antenna coords in lon/ lat in string format
    str_lon, str_lat, str_el = coords.ant_str

    # ITRF positions
    x, y, z = ant.positions.T.tolist()

    e, n, u = wgs84_to_enu(*coords.ant_rad, *coords.cofa)

    if args.u_file:
        u_file = args.u_file
    else:
        if ant.telescope == "MeerKAT":
            u_file = os.path.join(os.path.dirname(__file__),
                                  "data",
                                  "MeerKAT_WGS84_underlay.json")
        else:
            u_file = None

    # get underlay data
    if u_file:
        u_data = get_underlay_data(u_file)
        u_coords = get_antenna_coords(
            u_data.positions, u_data.telescope, offsets=u_data.offsets,
            cofa=u_data.cofa)

        re, rn, ru = wgs84_to_enu(*u_coords.ant_rad, *u_coords.cofa)
        r_names = u_data.names
    else:
        re, rn, ru = e, n, u
        r_names = ant.names

    source = dict(e=e, n=n, u=u, lon=str_lon, lat=str_lat, el=str_el,
                  x=x, y=y, z=z, indices=ant.indices, color=cmap,
                  name=ant.names, station=ant.stations, obs=ant.telescope,
                  re=re, rn=rn, ru=ru, rname=r_names,
                  rindices=np.where(np.isin(r_names, ant.names))[0])

    plot = plot_antennas(source, ms_name)

    save(plot)

    print(f"Antenna plot at: {o_file}")
