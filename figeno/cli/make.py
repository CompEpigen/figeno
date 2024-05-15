from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import traceback
from figeno import figeno_make
from figeno.utils import KnownException

def main(args):
    try:
        warnings=[]
        figeno_make(config_file=args.config,warnings=warnings)
        warning="\n".join(warnings)
        if warning!="":
            print("The figure was successfully generated, but with the following warnings:")
            print(warning)
            
    except KnownException as e:
        print("An error occured:")
        print(str(e))
    except Exception as e:
        print("An error occured.")
        print(traceback.format_exc())



def argparser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,add_help=False)
    parser.add_argument("config",type=str, help="Path to the config file (json).")
    return parser