from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(args):
    if not args.webview:
        from figeno.gui import gui_browser
        gui_browser.main(args)
    else:
        from figeno.gui import gui_webview
        gui_webview.main(args)

def argparser():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,add_help=False)
    parser.add_argument('-w','--webview', dest="webview",action='store_true',help="If set, start the GUI in webview (using pywebview) instead of in the browser..")
    parser.add_argument("-p","--port",type=int,default=5000, help="Port, only used in browser mode.")
    parser.add_argument('--debug', action='store_true',help="If set, will provide more debugging information.")
    parser.set_defaults(debug=False)
    parser.set_defaults(webview=False)
    return parser