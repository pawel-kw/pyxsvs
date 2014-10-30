Python XSVS code
================

Initialized during the sc3776 beamtime at ID10 (ESRF). Currently rewritten.
The basic principles of Speckle Visibility Spectroscopy are described in Bandyopadhyay et. al [1]

[1] Bandyopadhyay, R. et. al, Rev. Sci. Instr. *76*, 093110 (2005)

Required packages
-----------------

 * [SciPy](http://www.scipy.org/) - for plotting and calculations.
 * [fabio](http://sourceforge.net/apps/trac/fable/wiki/fabio) - for reading detector 
    data files in different formats. The code is tested with v. 0.1.4.
 * [pyFAI](https://github.com/kif/pyFAI) - fast azimuthal integration library. The code is tested with version 0.10.2.
 * [lmfit](http://cars9.uchicago.edu/software/python/lmfit/) - flexible and easy 
    to use non-linear least square fitting utilities.

Optionally:

 * [ReportLab](http://www.reportlab.com/opensource/) - open-source engine for creating complex, data-driven PDF documents and custom vector graphics. Required for the reportResults script to work.


Installation
------------

In principle installation should be as easy as typing:

    $ python setup.py install

Python setuptools should take care of the missing dependencies. In practice it's better
to install them by hand, epsecially pyFAI. 

Usage
-----

The 'examples' folder containes an IPython notebook which briefly describes the usage of pyxsvs. It can be loaded by running IPython from within the directory:

    $ ipython notebook

If you don't want to play with the notebook, a PDF version of it is also available in the 'examples' folder.

Another way of testing the code is to use the provided example data set and input file. The data is stored under

    ./examples/data/xsvs-series/

The calculation can be run by executing:

    $ cd ./examples/analysis/xsvs/
    $ pyxsvs -i xsvs_input.txt

If ReportLab is installed on your system, you can then run 

    reportResults -i ./xsvs_input.txt

to see a nice PDF report containing all the produced figures with some basic description.
