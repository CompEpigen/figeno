from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import platform

def main(args):
    mode = args.mode
    if mode=="auto":
        if platform.system()=="Linux": mode="browser"
        else: mode="webview"

    if mode=="browser":
        from figeno.gui import gui_browser
        gui_browser.main(args)
    else:
        from figeno.gui import gui_webview
        gui_webview.main(args)

def argparser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,add_help=False)
    parser.add_argument("-m","--mode",type=str,default="auto", 
                        help="The GUI can either be viewed in a browser tab (--mode browser) or in a separate webview window (--mode webview).\
                            By default (--mode auto), will select browser for linux and webview for windows and mac.")
    parser.add_argument("-p","--port",type=int,default=5000, help="Port, only used in browser mode.")
    parser.add_argument('--debug', action='store_true',help="If set, will provide more debugging information.")
    parser.set_defaults(debug=False)
    return parser