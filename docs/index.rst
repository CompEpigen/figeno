
Welcome to figeno's documentation!
==================================



Figeno is a tool to create (epi)genomics figures.

Features
---------

* Large collection of tracks (bigwig, HiC, alignments, copy number, SV...)
* Graphical user interface
* Multi-region figures
* Highlight regions of interest
* Output figures in svg, pdf, eps or png

A figure is described by a JSON file, which can either be created directly or generated using figeno's graphical user interface. The figure can then be generated either from the GUI or from the command line with `figenomake config.json`.


Creating the config file
------------------------
The config file is made out of five sections: General, Output, Regions, Highlights and Tracks.




.. toctree::
   :maxdepth: 2
   
   content/installation
   content/usage
   content/describe_figure
   content/examples

