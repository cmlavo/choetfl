"""
Function definitions for model/data import, export, processing, reading and writing operations related to model setup.
"""

import cobra
from etfl.core import Enzyme
import io
import os
import pandas as pd
from pprint import pprint
from pytfa.io import load_thermoDB
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

    if verbose: print("")

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
        pprint(model_validation)
    
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

    # read the .csv file, with advanced parsing so that the delimiter is ignored if surrounded by quotes
    df = pd.read_csv(io.StringIO(content),header=0, delimiter=',', quotechar="'")
    
    # set column names
    df.columns = ['prot_id','gene_id','component_id','product_id']
    
    # parse each row from text delimitation format to list or dict format based on the column
    df['gene_id'] = df['gene_id'].apply(lambda x: utils.parse_set(x)) # produces list from ';' delimiter
    df['component_id'] = df['component_id'].apply(lambda x: utils.parse_dict(x)) # produces dict from ';' delimiter and ':' key-value pair
    df['product_id'] = df['product_id'].apply(lambda x: utils.parse_set(x)) # produces list from ';' delimiter

    if verbose: print(f"{file_name} successfully imported [{df.shape[0]} rows x {df.shape[1]} columns].")

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

    if verbose: print(f"Merging {len(dataframes)} protein complex dataframes...")
    
    # Merge all dataframes into one
    merged_df = pd.concat(dataframes, ignore_index=True)
    initial_rows = merged_df.shape[0]

    # If items from the first row appear more than once, take the contents of the other columns (OrderedSets or dicts) and merge them while removing duplicates
    merged_df = merged_df.groupby('prot_id').agg({
        'gene_id': lambda x: utils.merge_set(x),
        'component_id': lambda x: utils.merge_dict(x, multiples_treatment),
        'product_id': lambda x: utils.merge_set(x)
    }).reset_index()

    final_rows = merged_df.shape[0]
    if verbose: 
        print(f"Merging operation removed {initial_rows-final_rows} rows.")
        print(f"Final merged protein complex dataframe: [{final_rows} rows x {merged_df.shape[1]} columns].")

    return merged_df

def assoc_prot_enz(dataframe: pd.DataFrame, k_constants: dict, verbose: bool = False) -> dict :
    """
    Create a dictionary associating each enzyme ID with its corresponding etfl.core.Enzyme object.
    Args:
        dataframe (pd.DataFrame): Dataframe containing enzyme data.
        k_constants (dict): Dictionary containing the k values configurations.
        verbose (bool): If True, prints status messages. Default is False.
    Returns:
        enzyme_dict (dict): Dictionary with protein complex IDs as keys and etfl.core.Enzyme objects as values.
    Raises:
        KeyError: If there is a key error in the product-gene or gene-component association.
    """

    enzyme_dict = {}
    count = 0

    for index, row in dataframe.iterrows():
        prot_id = row['prot_id']            # enzyme ID as a string
        gene_id = row['gene_id']            # OrderedSet of gene IDs
        component_id = row['component_id']  # dict of component IDs (with suffix) and their coefficients
        product_id = row['product_id']      # OrderedSet of product IDs (contains suffix)
        
        # for given protein compex, create a dictionary relating products to genes (by order of appearance). Discards any mismatches in numbers of products and genes   
        product_gene_assoc_dict = utils.assoc_gene_product(prot_id, product_id, gene_id, verbose) # for given protein complex, dictionary of product:gene pairs

        # create a dictionary for the composition of the protein complex which relates the gene ID for a given component to the component's coefficient
        composition = {product_gene_assoc_dict[product]: coeff for product, coeff in component_id.items() if product in product_gene_assoc_dict} # for all genes associated to a protein, dictionary of gene:component pairs
        
        if composition == {}:
            if verbose: print(f"Warning: No valid composition found for enzyme '{prot_id}'. Check for missing entries.")
            count += 1
            continue # do not create an Enzyme object for this enzyme if there is no valid composition

        enzyme = Enzyme(id = prot_id, kcat = k_constants.k_cat_default, kdeg = k_constants.k_deg_enz, composition = composition)
        enzyme_dict[prot_id] = enzyme

    if verbose: print(f"{len(enzyme_dict)} enzymes successfully associated with protein complexes. {count} enzymes discarded due to missing compositions.")

    return enzyme_dict

def import_kcat_data(file_name: str, default_kcat: float, verbose: bool = False) -> pd.DataFrame:
    """
    Imports kcat data from a .csv file with columns [Enzyme ID, kcat, kcat_bwd (sparse)] into a dataframe. Sets any
    missing k_cat_fwd values as the default value, and any missing k_cat_bwd values as the k_cat_fwd value.
    Args:
        file_name (str): The name of the .csv file to be imported.
        verbose (bool): If True, prints the import status. Default is False.
    Returns:    
        df (pd.DataFrame): A dataframe containing the kcat data.
    """

    file_path = os.path.join(os.getcwd(), "data/", file_name)

    assert os.path.exists(file_path), f"K_cat data file {file_name} does not exist in the data/ directory."
    assert file_path.endswith('.csv'), f"K_cat data file {file_name} must be a .csv file."

    df = pd.read_csv(file_path, header=0, delimiter=',', quotechar="'")

    # set column names and ensure k values are numerical
    df.columns = ['prot_id','k_cat_fwd','k_cat_bwd']

    df['k_cat_fwd'] = pd.to_numeric(df['k_cat_fwd'], errors='coerce').astype(float)
    df['k_cat_bwd'] = pd.to_numeric(df['k_cat_bwd'], errors='coerce').astype(float)

    # set any missing k_cat values to the default value
    df['k_cat_fwd'] = df['k_cat_fwd'].fillna(default_kcat)
    # set any missing k_cat_bwd values to the k_cat value
    df['k_cat_bwd'] = df['k_cat_bwd'].fillna(df['k_cat_fwd'])

    if verbose: print(f"{file_name} successfully imported [{df.shape[0]} rows x {df.shape[1]} columns].")

    return df

def import_thermo_data(file_name: str, verbose: bool = False) -> dict:
    """
    Imports thermodynamic data from a .thermodb file. The file must be located in the data/ directory.
    """

    file_path = os.path.join(os.getcwd(), "data/", file_name)
    assert os.path.exists(file_path), f"Thermodynamic data file {file_name} does not exist in the data/ directory."    

    thermo_db = load_thermoDB(file_path, verbose=verbose)

    if verbose: print(f"Successfully imported thermodynamic data from {file_name}.")

    return thermo_db