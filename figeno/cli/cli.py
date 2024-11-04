from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from figeno.cli import gui, init,make

__version__ = "1.6.1"

def main():
    parser = ArgumentParser("figeno",formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v','--version',action='version',version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(title="subcommands",description="valid commands",help="additional help",dest="command")
    subparsers.required = True

    for module in ["init","gui","make"]:
        mod = globals()[module]
        p = subparsers.add_parser(module,parents=[mod.argparser()])
        p.set_defaults(func=mod.main)

    args = parser.parse_args()
    args.func(args)