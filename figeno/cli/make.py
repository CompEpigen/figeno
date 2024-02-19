from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from figeno import figeno_make

def main(args):
    figeno_make(config_file=args.c)


def argparser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,add_help=False)
    parser.add_argument("-c","--config",type=str, help="Path to the config file (json).")
    parser.set_defaults(debug=False)
    return parser