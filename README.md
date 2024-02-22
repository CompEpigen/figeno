# figeno

Tool for making genomics figures.

## Features:
- Large collection of tracks (bigwig, HiC, alignments, copy number, SV…)
- Graphical user interface
- Multi-region figures with interactions across regions
- Highlight regions of interest
- Output figures in svg, pdf, eps or png

## Quick start
### Linux, MacOS
In an environment with python>=3.7:
```
pip install figeno
figeno gui
```
This will install figeno and run the graphical user interface (GUI). From the GUI, you can configure the figure and generate it, as well as save the JSON config file which fully defines the figure. The GUI is optional and you can instead use the command line interface: use [figeno init](https://figeno.readthedocs.io/en/latest/content/usage.html#figeno-init) to initialize a config file, edit the config file manually, and generate the figure with [figeno make](https://figeno.readthedocs.io/en/latest/content/usage.html#figeno-make).

### Windows
Download figeno_windows.zip from https://github.com/CompEpigen/figeno/releases/latest, unzip it and launch the graphical user interface by executing `figeno.exe`.

## Documentation
For more information on how to use figeno, please read the documentation at: 
https://figeno.readthedocs.io/en/latest/
  
