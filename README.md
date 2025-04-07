# choETFL
Applying the ETFL framework to the CHO genome-scale metabolic model.

## Setup Guide

This code requires [pyTFA](https://github.com/EPFL-LCSB/pytfa) and [ETFL](https://github.com/EPFL-LCSB/etfl) to function properly. To begin, clone a fork of these repos as well as a fork of this repo to a local folder with the following structure:

    └───src_folder
        └───pytfa
        └───etfl
        └───choetfl

It is recommended that you not clone the main branch of your personal forks, but to create and clone an unstable development branch of these forks. You can then preserve the main branch of your fork as a synchronisation branch to keep up with any updates made to these packages.

Next, navigate to ``choetfl`` and create the conda virtual environment using the ``requirements.yml`` file provided with this code. This step requires having Anaconda installed on your machine, see [installation instructions](https://www.anaconda.com/docs/getting-started/anaconda/install). From the ``src_folder``:

    cd choetfl
    conda env create -f requirements.yml
    conda activate choetfl_env

Next, navigate to pyTFA and install it as a package. Repeat this step for ETFL.

    cd ../pytfa
    pip install -e .
    cd ../etfl
    pip install -e .

Finally, an MILP solver such as Gurobi or CPLEX is required. Follow relevant installation instructions for the Python installation of your solver of choice. For example, if using Gurobi:

    conda install gurobi::gurobi

### Verifying Installation

You can run the ``test_installation.py`` script to ensure all libraries have successfully been installed with the correct version. You should expect the following output:

    Testing library versions:
    Python - 3.8.20
    NumPy - 1.19.5
    Pandas - 1.4.4
    Cobra - 0.24.0
    PyTFA - 0.9.4
    ETFL - 0.0.2
    GurobiPy - 11.0.3 (This will vary depending on the installed solver)

## Running the Code

