"""
This script executes all steps for the creation of the choETFL model.
"""

import os
import sys
from pathlib import Path

# Adding the parent directory to the system path
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))

from contextlib import redirect_stdout, redirect_stderr

import model_io
import model_ops
import utils

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
    baseline_model, _ = model_io.import_xml_model(config.baseline.file_name, config.solver, config.verbose)

    # Sanitizing metabolite and reaction names so that any IDs starting with a digit are prefixed with an underscore.
    #sanitized_model = model_ops.sanitize_model(baseline_model, config.verbose) TODO: Check if this is necessary
    with open(os.devnull, 'w') as fnull:
        with redirect_stdout(fnull), redirect_stderr(fnull):            
            sanitized_model = baseline_model.copy() # TODO: remove if unnecessary, replace sanitized_model with baseline_model in subsequent steps

    # Extracting the mass ratios of each biomass building block (BBB) from the biomass reaction of the baseline model
    mass_ratios = model_ops.compute_BBB_mass_ratios(sanitized_model, config.baseline.biomass_rxn_id, config.baseline.BBBs_params, config.verbose)
    
    ### MODEL MODIFICATION ###
    
    # Removing growth associated maintenance (GAM) cost of ATP from the biomass reaction; cost is later added explicitly
    GAM_modified_model = model_ops.modify_GAM(sanitized_model, config.baseline.biomass_rxn_id, config.baseline.BBBs_params, config.etfl.GAM_params, config.verbose)

    # Developing coupling dictionary of (Reaction:Enzyme list) pairs
    coupling_dict = model_ops.develop_coupling_dict(GAM_modified_model, config.etfl.enz_coupl_params, config.verbose)

if __name__ == "__main__":
    main("testing.yml")