from argparse import ArgumentParser
from Owlcat import PlotAntennas as pants


def get_arguments():
    """Command line arguments"""

    parser = ArgumentParser(usage="""%(prog)s: [options] """,
                            description="Plot antenna positions")

    parser.add_argument("--ms", dest="ms_name", required=True, metavar='',
                        help="MS for antennas. This is required",
                        default=None)
    parser.add_argument("-c", "--colour", dest="cmap",  metavar='', type=str,
                        help="""Use some color for antenna outline. Can be 
                        specified as colour name or hex / rgb / rgba. If
                        specified as hex value or rgb, use quotes. 
                        e.g'#1ECBE1' or 'rgb(255, 0, 0)'. """, default=None)
    parser.add_argument("--cofa", dest="cofa", metavar='',
                        type=str,
                        help="""Coordinates for the centrE of array in the
                        MS. Should be specified as a comma separated string
                        in the order lon,lat,elevation. Lon and Lat values
                        should be in degrees, while elevation should be in
                        metres. If no elevation is available, set the value
                        of 0.""")
    parser.add_argument("-o", "--output", dest="o_file",  metavar='',
                        type=str,
                        help="""Name to give output html file including the
                        file extension. e.g. 'test.html'. Default is
                        msName.html""", default=None)
    parser.add_argument("-uf", "--underlay-file", dest="u_file", metavar='',
                        type=str,
                        help="File containing entire array layout coordinates.",
                        default=None)

    return parser


if __name__ == "__main__":

    parser = get_arguments()
    args = parser.parse_args()
    pants.main(args)
