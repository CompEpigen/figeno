
Usage
==================================

Basics
-------

Figeno generates a figure based on a config file in JSON format. The config file describes where the output file will be stored, which regions will be shown, which tracks will be displayed... (see :ref:`Describing the figure`). 

The easiest way to create the config file is to use the graphical user interface, either by running ``figeno gui`` if you installed figeno as a python package, or by executing the executable if you downloaded the binaries. From the GUI, you can directly generate the figure, but also save the config file, in case you want to remake the figures later with slightly different parameters.

Alternatively, you can write the config file manually in a text editor, and then run ``figeno make -c /path/to/config.json`` to generate the figure. You can also initialize a config file which already contains some tracks using ``figeno init``.

    
Command line interace (CLI)
---------------------------

figeno init 
^^^^^^^^^^^

Initialize a config file, which will then have to be completed manually. Predefined templates are provided, or you can also provide a list of tracks and regions.

.. code:: bash

  figeno init -o config.json [--template TEMPLATE] [--tracks TRACKS] [--regions REGIONS] [--highlights HIGHLIGHTS]
  
Parameters:

* ``-o``, ``--output``: Path where the configuration file will be stored.
* ``--template``: Used to initialize the track with preset tracks and some options (see :ref:`Templates`). Possible values:

  * bigwig. Will have the following tracks: bigwig, genes, chr_axis.
  * hic. Will have the following tracks: hic, genes, chr_axis.
  * asm. Will have the following tracks: alignments (with split by haplotype and color by basemod), basemod_freq, genes, chr_axis.
  * wgs_chr. Will have the following tracks: sv, copynumber, chr_axis. Will also set margin_above to 0.
  * wgs_circos. Will have the following tracks: sv, copynumber, chr_axis. Will also set margin_above to 0, layout to circular and set the regions to all chromosomes.
  
* ``--tracks``: comma-separated list of tracks, eg: bigwig,genes,chr_axis. 

* ``--regions``: comma-separated list of regions, eg: 7:156790000-156820000,12:11900000-12100000

* ``--highlights``: comma-separated list of highlights (same format as regions).
  
This is a normal text paragraph. The next paragraph is a code sample::

   It is not processed in any way, except
   that the indentation is removed.

   It can span multiple lines.

This is a normal text paragraph again.
   
figeno make
^^^^^^^^^^^

Generate a figure described by a config file.

.. code:: bash

  figeno make -c config.json
  
Parameters:

* ``-c``, ``--config``: config file.


figeno gui
^^^^^^^^^^^

Start the graphical user interface.
   
.. code:: bash

  figeno gui [-p PORT] [--debug]
  
Parameters:

* ``-p``,``--port``. Port for the local server (default: 5000).

* ``--debug``: if set, will print more information to the terminal.
  

  

   





