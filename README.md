# figeno

Tool for making genomics figures.

![figeno](docs/content/images/figeno.png)
Example figures generated with figeno. Left: allele-specific methylation with nanopore data. Right: HiC data across a structural rearrangement.

## Features
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

## Test data
Example input files to test figeno are provided in [test_data](https://github.com/CompEpigen/figeno/tree/main/test_data).

## Bugs and improvements
If you encounter a bug or would like to have a new feature added, please do not hesitate to [raise an issue](https://github.com/CompEpigen/figeno/issues/new) or to [contact me directly](https://www.dkfz.de/en/CanEpi/staff/kontakt/Sollier_Etienne.php).
  
