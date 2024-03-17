
Usage
==================================

Basics
-------

Figeno generates a figure based on a config file in JSON format. The config file describes where the output file will be stored, which regions will be shown, which tracks will be displayed... (see :ref:`Describing the figure` and `examples of config files <https://github.com/CompEpigen/figeno/tree/main/test_data>`_). 

The easiest way to create the config file is to use the graphical user interface (GUI), either by running ``figeno gui`` if you installed figeno as a python package, or by executing the executable if you downloaded the binaries. From the GUI, you can directly generate the figure, but also save the config file, in case you want to remake the figures later with slightly different parameters.

Alternatively, you can write the config file manually in a text editor, and then run ``figeno make /path/to/config.json`` to generate the figure. You can also initialize a config file which already contains some tracks using ``figeno init``.

    
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
  * wgs_circos. Will have the following tracks: sv, copynumber, chr_axis. Will also set margin_above to 0, max_cn to 3.9 and disable the vertical grid for the copynumber track, set layout to circular and set the regions to all chromosomes.
  
* ``--tracks``: comma-separated list of tracks, eg: bigwig,genes,chr_axis. 

* ``--regions``: comma-separated list of regions, eg: 7:156790000-156820000,12:11900000-12100000

* ``--highlights``: comma-separated list of highlights (same format as regions).
  
   
figeno make
^^^^^^^^^^^

Generate a figure described by a config file.

.. code:: bash

  figeno make config.json
  
Parameters:

* config file in json format.


figeno gui
^^^^^^^^^^^

Start the graphical user interface.
   
.. code:: bash

  figeno gui [--webview] [-p PORT] [--debug]
  
Parameters:

* ``-w``, ``--webview``. If set, will use pywebview to render the GUI. Otherwise, the GUI can be viewed in the browser (at localhost:5000 by default).

* ``-p``, ``--port``. Port for the local server, in case the browser mode is used (default: 5000).

* ``--debug``: If set, will print more information to the terminal.

.. warning::
  The webview mode works for linux, windows and mac, but for linux you will need to install additional dependencies (see https://pywebview.flowrl.com/guide/installation.html#linux).
  

Graphical user interface (GUI)
------------------------------

The GUI can be started with ``figeno gui`` from the command line, or by launching the executable for windows. It can be used to easily edit a config file. Required parameters which have not been filled in yet are highlighted in orange. Once you have finished describing the figure, you can click on "Generate figure" to generate it. You can also save the config file to a json file by clicking on "Save config" if you want to edit it later, in which case you can then load it again with "Load config". You can also combine the CLI and the GUI, for example by creating a config file with the GUI, saving it, and then using ``figeno make`` to generate the figure.

Python API
-----------

You can also import figeno as a python module, and give ``figeno_make`` the config file as a python dictionary.


.. code:: python

  import figeno_make from figeno
  
  config={"general":{"reference":"hg19","layout":"horizontal"}}
  config["output"] = {"file":"figure.svg","dpi":200,"width":180}
  config["regions"] = [{"chr":"17","start":7000000,"end":7500000}]
  config["tracks"] = [{"type":"genes"}, {"type":"chr_axis"}]
  figeno_make(config)
  

   





