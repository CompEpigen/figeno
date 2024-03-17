
Installation
==================================

    
Linux, MacOS
^^^^^^^^^^^^

In an environment with python>=3.7:

.. code:: bash
	
	pip install figeno


Windows
^^^^^^^

Figeno has some dependencies which do not officially support Windows (pysam, pybigwig), so it is difficult to install as a python package on Windows. To make it easier for windows users to use figeno, we provide windows binaries for figeno. Download figeno_windows.zip from https://github.com/CompEpigen/figeno/releases/latest, extract it, and start the figeno GUI by double-clicking on figene.exe.


Building from source
^^^^^^^^^^^^^^^^^^^^

It is easier to install figeno from pip or with the windows binaries. If you want to build figeno from the source code and want to use the GUI, you will first need to build the react app, which requires `nodejs <https://nodejs.org/en>`_ to be installed (giving you access to the npm command), and then install the python package.

.. code:: bash

	cd figeno/gui
	npm install
	npm run build
	cd ../..
	pip install .

