"""
This script executes all steps for the creation of the choETFL model.
"""

import os
import sys

# Adding the parent directory to the system path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'project_code'))
sys.path.append(parent_dir)

import model_io
import model_ops
from project_code import utils

def main(config: str = "default.yml"):
    """
    Main method for setting up the choETFL model.
    Args:
        config (str): Name of the configuration file.
    """

    ### CONFIGURATION ###

    # Loading the configuration file
    config = utils.load_config(config)

    ### BASELINE MODEL ###

    # Importing the baseline iCHO2291 model in .xml format from which choETFL is constructed.
    baseline_model, model_validation = model_io.import_json_model(config.baseline.file_name, config.solver, config.baseline.verbose)

    # Sanitizing metabolite and reaction names
    #baseline_model = model_ops.sanitize_model(baseline_model)

    ### 
    

if __name__ == "__main__":
    main("testing.yml")