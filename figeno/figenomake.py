from figeno.tracks_plot import tracks_plot
import sys


def main():
    if len(sys.argv)<2: raise Exception("Must provide the config file.")
    tp = tracks_plot(config_file=sys.argv[1])
    tp.draw()

if __name__=="__main__":
    main()
    