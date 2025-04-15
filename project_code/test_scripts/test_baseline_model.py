"""
Script which executes various tests to ensure the integrity and functionality of the baseline iCHO2291 model.
"""

# Import baseline model
# Run checks on biomass reaction, subreactions
#   - Check that all molar weights are available, if not then suggest edits in config file
# 

import os
import sys
from pathlib import Path

# Adding the parent directory to the system path
parent_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(parent_dir))

from model_setup import model_io
from model_setup import model_ops
import utils

config = utils.load_config("testing.yml")
model, _ = model_io.import_xml_model(config.baseline.file_name, config.solver, config.verbose)

reactions_view = [reaction for reaction in model.reactions if model.metabolites.get_by_id('adp[c]') in reaction.reactants]

breakpoint()