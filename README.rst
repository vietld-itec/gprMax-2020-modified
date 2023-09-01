.. image:: https://readthedocs.org/projects/gprmax/badge/?version=latest
    :target: http://docs.gprmax.com/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: docs/source/images/gprMax_logo_small.png
    :target: http://www.gprmax.com

***************
Getting Started
***************

Installation
============

1. Install Python, required Python packages, and get gprMax source
------------------------------------------------------------------
* `Download Anaconda Python 3.8 <https://repo.anaconda.com/archive/Anaconda3-2020.11-Windows-x86_64.exe>`_.
* Open a Terminal Command Prompt (Windows) and run the following commands:

.. code-block:: bash

    $ git clone https://github.com/vietld-itec/gprMax-rebuild
    $ cd gprMax-rebuild
    $ conda env create -f conda_env.yml

This will make sure conda is up-to-date, install Git, get the latest gprMax source code from GitHub, and create an environment for gprMax with all the necessary Python packages.


2. Install a C compiler which supports OpenMP
---------------------------------------------

* Download and install `Microsoft Visual C++ 2015 Build Tools <BuildTools_Full.exe>`_.

3. Build and install gprMax
---------------------------

Once you have installed the aforementioned tools follow these steps to build and install gprMax:

* Open a Terminal Command Prompt (Windows), navigate into the top-level gprMax directory, and if it is not already active, activate the gprMax conda environment :code:`conda activate gprMax`. Run the following commands:

.. code-block:: bash

    (gprMax)$ python setup.py build
    (gprMax)$ python setup.py install

**You are now ready to proceed to running gprMax.**

If you have problems with building gprMax on Microsoft Windows, you may need to add :code:`C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin` to your path environment variable.

Running gprMax
==============

gprMax is designed as a Python package, i.e. a namespace which can contain multiple packages and modules, much like a directory.

Open a Terminal (Linux/macOS) or Command Prompt (Windows), navigate into the top-level gprMax directory, and if it is not already active, activate the gprMax conda environment :code:`conda activate gprMax`.

Basic usage of gprMax is:

.. code-block:: bash

    (gprMax)$ python -m gprMax path_to/name_of_input_file

For example to run one of the test models:

.. code-block:: bash

    (gprMax)$ python -m gprMax user_models/cylinder_Ascan_2D.in

When the simulation is complete you can plot the A-scan using:

.. code-block:: bash

    (gprMax)$ python -m tools.plot_Ascan user_models/cylinder_Ascan_2D.out

Your results should like those from the A-scan from the metal cylinder example in `introductory/basic 2D models section <http://docs.gprmax.com/en/latest/examples_simple_2D.html#view-the-results>`_

When you are finished using gprMax, the conda environment can be deactivated using :code:`conda deactivate`.

Optional command line arguments
-------------------------------

====================== ========= ===========
Argument name          Type      Description
====================== ========= ===========
``-n``                 integer   number of times to run the input file. This option can be used to run a series of models, e.g. to create a B-scan with 60 traces: ``(gprMax)$ python -m gprMax user_models/cylinder_Bscan_2D.in -n 60``
``-gpu``               flag/list flag to use NVIDIA GPU or list of NVIDIA GPU device ID(s) for specific GPU card(s), e.g. ``-gpu 0 1``
``-restart``           integer   model number to start/restart simulation from. It would typically be used to restart a series of models from a specific model number, with the ``-n`` argument, e.g. to restart from A-scan 45 when creating a B-scan with 60 traces: ``(gprMax)$ python -m gprMax user_models/cylinder_Bscan_2D.in -n 15 -restart 45``
``-task``              integer   task identifier (model number) when running simulation as a job array on `Open Grid Scheduler/Grid Engine <http://gridscheduler.sourceforge.net/index.html>`_. For further details see the `parallel performance section of the User Guide <http://docs.gprmax.com/en/latest/openmp_mpi.html>`_
``-mpi``               integer   number of Message Passing Interface (MPI) tasks, i.e. master + workers, for MPI task farm. This option is most usefully combined with ``-n`` to allow individual models to be farmed out using a MPI task farm, e.g. to create a B-scan with 60 traces and use MPI to farm out each trace: ``(gprMax)$ python -m gprMax user_models/cylinder_Bscan_2D.in -n 60 -mpi 61``. For further details see the `parallel performance section of the User Guide <http://docs.gprmax.com/en/latest/openmp_mpi.html>`_
``--mpi-no-spawn``     flag      use MPI task farm without spawn mechanism. For further details see the `parallel performance section of the User Guide <http://docs.gprmax.com/en/latest/openmp_mpi.html>`_
``-benchmark``         flag      switch on benchmarking mode. This can be used to benchmark the threading (parallel) performance of gprMax on different hardware. For further details see the `benchmarking section of the User Guide <http://docs.gprmax.com/en/latest/benchmarking.html>`_
``--geometry-only``    flag      build a model and produce any geometry views but do not run the simulation, e.g. to check the geometry of a model is correct: ``(gprMax)$ python -m gprMax user_models/heterogeneous_soil.in --geometry-only``
``--geometry-fixed``   flag      run a series of models where the geometry does not change between models, e.g. a B-scan where *only* the position of simple sources and receivers, moved using ``#src_steps`` and ``#rx_steps``, changes between models.
``--opt-taguchi``      flag      run a series of models using an optimisation process based on Taguchi's method. For further details see the `user libraries section of the User Guide <http://docs.gprmax.com/en/latest/user_libs_opt_taguchi.html>`_
``--write-processed``  flag      write another input file after any Python code and include commands in the original input file have been processed. Useful for checking that any Python code is being correctly processed into gprMax commands.
``-h`` or ``--help``   flag      used to get help on command line options.
====================== ========= ===========

Updating gprMax
===============

* Open a Terminal (Linux/macOS) or Command Prompt (Windows), navigate into the top-level gprMax directory, and if it is not already active, activate the gprMax conda environment :code:`conda activate gprMax`. Run the following commands:

.. code-block:: bash

    (gprMax)$ git pull
    (gprMax)$ python setup.py cleanall
    (gprMax)$ python setup.py build
    (gprMax)$ python setup.py install

This will pull the most recent gprMax source code form GitHub, remove/clean previously built modules, and then build and install the latest version of gprMax.


Updating conda and Python packages
----------------------------------

Periodically you should update conda and the required Python packages. With the gprMax environment deactivated and from the top-level gprMax directory, run the following commands:

.. code-block:: bash

    $ conda update conda
    $ conda env update -f conda_env.yml


How to use plot_A_scan_raw.py and plot_Bscan_gain.py (Developed by Viet Le @ 2020)
===============
1. plot_A_scan_raw.py
------------------------------------------------------------------
a/	Export raw data and normalized data ( *.csv files): Plot Ascan, export rawdata and normalized data ( 2 csv files)

.. code-block:: bash

    python –m tools.plot_Ascan_raw outputfile.out –outputs Ez –rawdata
    
b/	Plotting Ascan from rawdata. I won’t export data (not thing, don’t have *.csv file, this feature like plot_Ascan.py from gprMax package)

.. code-block:: bash

    python –m tools.plot_Ascan_raw outputfile.out –outputs Ez
    
c/ Plotting Ascan from normalized data

.. code-block:: bash

    python –m tools.plot_Ascan_raw outputfile.out –outputs Ez -plotnorm
    
2. plot_B_scan_gain.py
------------------------------------------------------------------
a/	Plotting Bscan with raw data (time domain)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out
    
b/	Plotting Bscan with normalized data (time domain)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out -norm
    
c/	Plotting Bscan with er (equivalent relative dielectric constant) with raw data (including depth axis)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out –er 4.8
    
d/	Plotting Bscan with er (equivalent relative dielectric constant) with normalized data (including depth axis)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out –er 4.8 -norm
    
e/	Plotting Bscan with er (equivalent relative dielectric constant) with raw data and also apply gain function (including depth axis)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out –er 4.8 –gmin 1 –gmax 50
    
f/	Plotting Bscan with er (equivalent relative dielectric constant) with normalized data and also apply gain function (including depth axis)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out –er 4.8 -norm –gmin 1 –gmax 50
    
g/	Plotting Bscan with raw data and also apply gain function (time domain)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out –gmin 1 –gmax 50
    
h/	Plotting Bscan with normalized data and also apply gain function (time domain)

.. code-block:: bash

    python –m tools.plot_Bscan_gain Ez outputfile.out –gmin 1 –gmax 50 -norm
    
(Contact to Author: viet.xd.bkdn@gmail.com / Kakaotalk ID: vietld1991)
