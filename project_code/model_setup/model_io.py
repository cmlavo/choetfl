"""
Function definitions for import, export, reading and writing operations related to model setup.
"""

import cobra
import os
import pprint
import time
from typing import Tuple

def import_xml_model(model_name: str, solver: str = 'gurobi', verbose: bool = False) -> Tuple[cobra.Model, str]:
    """
    Imports a .xml model located in the /models directory based on the file name provided and the selected solver.

    Args:
        model_name (str): The name of the .json model file to be imported.
        solver (str): The solver to be used. Default is 'gurobi'.
        verbose (bool): If True, prints the validation status of the model. Default is False.
    Returns:
        cobra_model (cobra.Model): The imported COBRA model.
        model_validation (str): The validation status of the model.
    """

    model_dir = os.path.join(os.getcwd(), "models/", model_name)

    assert os.path.exists(model_dir), f"Model file {model_name} does not exist in the /models directory."
    assert model_name.endswith('.xml'), "Model must be an .xml file"

    if verbose: print(f"Importing {model_name}...")

    start_time = time.perf_counter()
    cobra_model, model_validation = cobra.io.validate_sbml_model(model_dir)
    end_time = time.perf_counter()
    execution_time = end_time - start_time
    if verbose: print(f"{model_name} import complete ({execution_time:.2f} seconds).")

    cobra_model.solver = solver    

    if verbose: 
        print("Errors:")
        pprint.pprint(model_validation)
    
    return cobra_model, model_validation