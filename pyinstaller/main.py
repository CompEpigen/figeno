import argparse
import figeno
import figeno.cli.gui

args=argparse.ArgumentParser().parse_args()
args.webview=True
args.port=5000
args.debug=False
figeno.cli.gui.main(args)

