from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from figeno.gui import gui

def main(args):
    gui.main(args)


def argparser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,add_help=False)
    parser.add_argument("-p","--port",type=int,default=5000, help="Port.")
    parser.add_argument('--debug', action='store_true')
    parser.set_defaults(debug=False)
    return parser