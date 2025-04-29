"""
Function definitions for general operations.
"""

from box import Box
import math
from ordered_set import OrderedSet
import os
import pandas as pd
from pprint import pprint
import yaml

def load_config(path: str = "default.yml") -> Box:
    """
    Load a config .yml file which specifies configurations for all scripts.
    Args:
        path (str): Path to the config .yml file. Default is "configs/default.yml".
    Return:
        config (box.Box): Box dictionary containing the configurations.
    """
    
    with open(os.path.join("configs/", path), "r") as f:
        config = yaml.safe_load(f)
    return Box(config)

def remove_unquoted_whitespace(content: str, quotechar: str = "'") -> str:
    """
    Remove all whitespace from a string, but only if not surrounded by quotes.
    Args:
        content (str): The string to be processed.
        quotechar (str): The character used for quoting. Default is '.
    Returns:
        str: The processed string with unquoted whitespace removed.
    """
    string_list = []
    in_quotes = False
    for char in content:
        if char == quotechar:
            in_quotes = not in_quotes
            string_list.append(char)
        elif char == " " and not in_quotes:
            continue
        else:
            string_list.append(char)
    return "".join(string_list)

def parse_set(content: str, delimiter: str = ';', quotechar: str = "'") -> OrderedSet:
    """
    Parse a string representation of a set, separated by a delimiter, into an actual set, ignoring occurences 
    of this delimiter within quotes. Does not remove whitespace, please ensure that this is done beforehand.
    Args:
        content (str): The string to be parsed. If Nan, returns None.
        delimiter (str): The character used to separate items in the set. Default is ';'.
        quotechar (str): The character used for quoting. Default is '.
    Returns:
        parsed_set (OrderedSet): The parsed set.
    """

    if isinstance(content, float) and math.isnan(content): return OrderedSet()

    parsed_set = OrderedSet()
    in_quote = False
    current_item = ""
    for char in content:        
        if char == quotechar:
            in_quote = not in_quote
        elif char == delimiter and not in_quote:
            parsed_set.add(current_item)
            current_item = ""
        else:
            current_item += char
    parsed_set.add(current_item)  # Add the last item

    return parsed_set

def parse_dict(content: str, delimiter: str = ';', associator: str = ':', quotechar: str = "'") -> list:
    """
    Parse a string representation of a dict of key (str) : value (float), separated by a delimiter and with key-value pairs 
    denoted by an associator, into an actual dict, ignoring occurences of this delimiter within quotes. Does not remove 
    whitespace, please ensure that this is done beforehand.
    Args:
        content (str): The string to be parsed. If Nan, returns None.
        delimiter (str): The character used to separate items in the dict. Default is ';'.
        associator (str): The character used to separate keys and values. Default is ':'.
        quotechar (str): The character used for quoting. Default is '.
    Returns:
        dict (dict): The parsed dict.
    """

    if isinstance(content, float) and math.isnan(content): return dict()

    parsed_dict = {}
    in_quote = False
    current_item = ""    
    for char in content:
        if char == quotechar:
            in_quote = not in_quote
        elif char == delimiter and not in_quote:
            if associator in current_item:
                key, value = current_item.split(associator, 1)
                parsed_dict[key] = float(value)
            current_item = ""
        else:
            current_item += char

    # Handle the last item after the loop
    if current_item:
        if associator in current_item:
            key, value = current_item.split(associator, 1)
            parsed_dict[key] = int(value)

    return parsed_dict

def merge_set(series: pd.Series) -> OrderedSet:
    """
    Merge a series of sets into a single set, removing duplicates.
    Args:
        series (pd.Series): The series to be merged.
    Returns:
        set (OrderedSet): The merged set.
    """
    merged_set = OrderedSet()
    for item in series:
        assert isinstance(item, OrderedSet), f"Item \'{item}\' is not an OrderedSet."
        merged_set.update(item)

    return merged_set

def merge_dict(series: pd.Series, multiples_treatment: int) -> dict:
    """
    Merge a series of dicts into a single dict, removing duplicate pairs. If a key appears more than once but with
    different values, values are treated in the specified way.
    Args:
        series (pd.Series): The series of dicts to be merged.
        multiples_treatment (int): The treatment for multiple values:
    0 - Keep the lowest value.  
    1 - Keep the highest value.  
    2 - Take the average value.  
    3 - Take the value which appears first in the series.  
    4 - Take the value which appears last in the series.  
    Returns:
        dict: The merged dict.
    """

    assert multiples_treatment in [0, 1, 2, 3, 4], f"Invalid multiples_treatment value: {multiples_treatment}. Must be an integer between 0 and 4."

    merged_dict = {}
    for item in series:
        assert isinstance(item, dict), f"Item \'{item}\' is not a dict."
        for key, value in item.items():
            if key in merged_dict and value != merged_dict[key]:
                print(f"Warning: Component \'{key}\' appears multiple times in the protein complex data with different coefficient values. Merging values using treatment method {multiples_treatment}.")
                if multiples_treatment   == 0: merged_dict[key] = min(merged_dict[key], value)
                elif multiples_treatment == 1: merged_dict[key] = max(merged_dict[key], value)
                elif multiples_treatment == 2: merged_dict[key] = (merged_dict[key] + value) / 2
                elif multiples_treatment == 3: pass
                elif multiples_treatment == 4: merged_dict[key] = value
            elif key in merged_dict and value == merged_dict[key]: pass
            else: merged_dict[key] = value

    return merged_dict

def assoc_gene_product(prot_id: str, product_id: OrderedSet, gene_id: OrderedSet, verbose: bool) -> dict:
    """
    Create a dictionary associating each product ID with its corresponding gene ID, by order of appearance
    in the OrderedSet. Discards any genes for which there is no product, or any product for which there is no gene.
    Args:
        prot_id (str): Protein complex ID as a string.
        product_id (OrderedSet): OrderedSet of product IDs.
        gene_id (OrderedSet): OrderedSet of gene IDs.
        verbose (bool): If True, prints warning messages. Default is False.
    Returns:
        dict: Dictionary with product IDs as keys and gene IDs as values.
    """

    # Check if each product ID has a corresponding gene ID; if not, warn about discarding the extra products or genes
    if verbose:
        if len(product_id) > len(gene_id):
            print(f"Warning: Number of products ({len(product_id)}) exceeds number of genes ({len(gene_id)}) for enzyme '{prot_id}'. Check for duplicate or missing entries.")
            print(f"Discarding products: {list(product_id[len(gene_id):])}")        
        elif len(product_id) < len(gene_id):
            print(f"Warning: Number of genes ({len(gene_id)}) exceeds number of products ({len(product_id)}) for enzyme '{prot_id}'. Check for duplicate or missing entries.")
            print(f"Discarding genes: {list(gene_id[len(product_id):])}")        

    # Return dictionary associating each product ID with its corresponding gene ID, discarding mismatches
    return dict(zip(product_id, gene_id))

def jaccard_score(set1: OrderedSet, set2: OrderedSet) -> float:
    """
    Calculate the Jaccard score between two sets.
    Args:
        set1 (OrderedSet): First set.
        set2 (OrderedSet): Second set.
    Returns:
        float: Jaccard score between the two sets.
    """

    intersection = set1 & set2
    union = set1 | set2
    if not union:
        return 1.0  # define Jaccard of two empty sets as 1
    return len(intersection) / len(union)

    