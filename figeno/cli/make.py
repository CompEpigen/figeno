from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from figeno import figeno_make

def main(args):
    figeno_make(config_file=args.config)


def argparser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,add_help=False)
    parser.add_argument("config",type=str, help="Path to the config file (json).")
    return parser