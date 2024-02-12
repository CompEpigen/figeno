.. figeno documentation master file, created by
   sphinx-quickstart on Thu Feb  8 16:49:35 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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

General
^^^^^^^

* **figure_layout**: 

  * **horizontal**: Draw all regions one after the other horizontally
  
  * **circular**: Draw all regions one after the other in a circle
  
  * **symmetric**: Draw the regions in two rows, such that the tracks are symmetric: the bottom row has the tracks in the normal order, but the top row has its tracks in the reverse order. This is mainly intended to show copy numbers and breakpoints, with sv as the topmost track.
  
  * **stacked**" Draw all regions horizontally, but instead of being next to each other horizontally, the different regions are stacked vertically.
* **reference**: Reference genome used. Files for hg19 and hg38 are provided. You can also choose a custom reference genome, but then you will need to provide the required files.

Output
^^^^^^

* **file**: Path to the output file. The format will be inferred from the file extension
* **dpi**: Higher values will result in higher resolution but larger file size. Even if the figure is saved as vector graphics, this parameter will still have an importance for the hic and alignments tracks because they are rasterized.
* **width**: Width of the figure in mm. The default is 180mm, which is a standard full-page figure.

Regions
^^^^^^^

Regions are defined by chr, start and end. If end<start, the region will be shown in reverse orientation (from right to left). The color attribute will only be used if a chr axis track with style "arrow" is used. At least one region must be provided. 



Highlights
^^^^^^^^^^

Optionally, highlight some areas defined by chr, start and end, for example peaks.

Tracks
^^^^^^


.. toctree::
   :maxdepth: 1
   
   content/tracks

