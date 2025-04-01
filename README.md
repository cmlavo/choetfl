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
    conda install -e .
    cd ../etfl
    conda install -e .

Finally, an MILP solver such as Gurobi or CPLEX is required. Follow relevant installation instructions for the Python installation of your solver of choice.

## Running the Code

