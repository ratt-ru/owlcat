from argparse import ArgumentParser
from Owlcat import PlotAntennas as pants


def get_arguments():
    """Command line arguments"""

    parser = ArgumentParser(usage="""%(prog)s: [options] """,
                            description="Plot antenna positions")

    parser.add_argument("-ms", dest="ms_name", required=True, metavar='',
                        help="MS for antennas. This is required",
                        default=None)
    parser.add_argument(
        "-cmap", dest="cmap",  metavar='',
        help="""Use some colormap for the antennas. The default colour map is
        Virids. Allowable colours can be found at:
        https://docs.bokeh.org/en/latest/docs/reference/palettes.html. 
        Note the case sensitivity.""",
        default=None)
    parser.add_argument("-o", dest="o_file",  metavar='',
                        help="""Name to give output html file including the file extension.
        e.g. 'test.html'. Default is msName.html""",
                        default=None)

    return parser


if __name__ == "__main__":

    parser = get_arguments()
    args = parser.parse_args()
    pants.main(args)
