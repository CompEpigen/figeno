
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

.. argparse::
   :module: figeno.cli.init
   :func: argparser
   :prog: figeno init
   
figeno make
^^^^^^^^^^^

.. argparse::
   :module: figeno.cli.make
   :func: argparser
   :prog: figeno make

figeno gui
^^^^^^^^^^^
   
.. argparse::
   :module: figeno.cli.gui
   :func: argparser
   :prog: figeno gui
   



Graphical user interface (GUI)
-------------------------------





