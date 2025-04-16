"""
Function definitions for model/data import, export, processing, reading and writing operations related to model setup.
"""

import cobra
import io
import os
import pandas as pd
import pprint
import time
from typing import Tuple
import utils

def import_xml_model(model_name: str, solver: str = 'gurobi', verbose: bool = False) -> Tuple[cobra.Model, str]:
    """
    Imports a .xml model located in the models/ directory based on the file name provided and the selected solver.

    Args:
        model_name (str): The name of the .xml model file to be imported.
        solver (str): The solver to be used. Default is 'gurobi'.
        verbose (bool): If True, prints the validation status of the model. Default is False.
    Returns:
        cobra_model (cobra.Model): The imported COBRA model.
        model_validation (str): The validation status of the model.
    """

    model_path = os.path.join(os.getcwd(), "models/", model_name)

    if verbose: print("\n")

    assert os.path.exists(model_path), f"Model file {model_name} does not exist in the models/ directory."
    assert model_name.endswith('.xml'), "Model must be an .xml file."

    if verbose: print(f"Importing {model_name}...")

    start_time = time.perf_counter()
    cobra_model, model_validation = cobra.io.validate_sbml_model(model_path)
    end_time = time.perf_counter()
    execution_time = end_time - start_time
    if verbose: print(f"{model_name} import complete ({execution_time:.2f} seconds).")

    cobra_model.solver = solver

    if verbose: 
        print("Errors:")
        pprint.pprint(model_validation)
    
    return cobra_model, model_validation

def export_model(model: cobra.Model, model_name: str, solver: str = 'gurobi', verbose: bool = False) -> None:
    """
    Exports a COBRA ETFL model to a .xml file.
    """
    pass

def import_protein_complex_data(file_name: str, verbose: bool = False) -> pd.DataFrame:
    """
    Imports protein complex data from a .csv file with columns [Product Name, Genes (list separated by ';'), 
    Component Coefficients (list separated by ';' of tuples associated by ':'), Gene Products (list separated by ';')] 
    separated by ','. Returns a dataframe with the data.
    Args:
        file_name (str): The name of the .csv file to be imported.
        verbose (bool): If True, prints the import status. Default is False.
    Returns:
        pd.DataFrame: A dataframe containing the protein complex data.
    """
    file_path = os.path.join(os.getcwd(), "data/", file_name)

    assert os.path.exists(file_path), f"Protein complex data file {file_name} does not exist in the data/ directory."
    assert file_path.endswith('.csv'), f"Protein complex data file {file_name} must be a .csv file."

    # change all quote characters from " to ' in the file, since .csv readers only accept one type of quote character
    # also strip all whitespace if not surrounded by quotes
    with open(file_path, 'r', encoding='utf-8') as f:
        content = utils.remove_unquoted_whitespace(f.read().replace('"', "'"))
        breakpoint()

    # read the .csv file, with advanced parsing so that the delimiter is ignored if surrounded by quotes
    df = pd.read_csv(io.StringIO(content),header=0, delimiter=',', quotechar="'")

    if verbose: print(f"{file_name} successfully imported [{df.shape[0]} rows x {df.shape[1]} columns].")
    
    # set column names
    df.columns = ['prot_id','gene_id','components_ids','products_ids']

    # parse each row from text delimitation format to list or dict format based on the column
    df['gene_id'] = df['gene_id'].apply(lambda x: utils.parse_set(x)) # produces list from ';' delimiter
    df['components_ids'] = df['components_ids'].apply(lambda x: utils.parse_dict(x)) # produces dict from ';' delimiter and ':' key-value pair
    df['products_ids'] = df['products_ids'].apply(lambda x: utils.parse_set(x)) # produces list from ';' delimiter

    return df

def merge_protein_complex_data(dataframes: list, multiples_treatment: int = 3, verbose: bool = False) -> pd.DataFrame:
    """
    Merges multiple protein complex dataframes loaded using import_protein_complex_data() into a single dataframe, 
    ensuring no duplicates exist. If a protein ID appears in multiple dataframes but the other columns differ, the 
    information is combined into one row.
    Args:
        dataframes (list): List of pd.Dataframe objects to be merged.
        multiples_treatment (int): Treatment for multiple values:
            0: Keep the lowest value
            1: Keep the highest value
            2: Keep the first value
            3: Keep the last value
        verbose (bool): If True, prints the merge status. Default is False.
    Returns:
        pd.DataFrame: Merged dataframe with unique, combined entries.
    """
    
    # Merge all dataframes into one
    merged_df = pd.concat(dataframes, ignore_index=True)

    # If items from the first row appear more than once, take the contents of the other columns (sets or dicts) and merge them while removing duplicates
    merged_df = merged_df.groupby('prot_id').agg({
        'gene_id': lambda x: utils.merge_set(x),
        'components_ids': lambda x: utils.merge_dict(x, multiples_treatment),
        'products_ids': lambda x: utils.merge_set(x)
    }).reset_index()

    return merged_df
