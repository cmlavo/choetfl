"""
Function definitions for model operations related to model setup.
"""

import cobra
import pprint

def sanitize_model(model: cobra.Model) -> cobra.Model:
    """
    Change all reaction and metabolite IDs that start with a digit to start with an underscore instead.
    Args:
        model (cobra.Model): The model to sanitize.
    Returns:
        sanitized_model (cobra.Model): The sanitized model.
    """

    for met in model.metabolites:
        if met.id[0].isdigit():
            met.id = '_' + met.id
    for rxn in model.reactions:
        if rxn.id[0].isdigit():
            rxn.id = '_' + rxn.id
    model.repair()

    return model

def get_MWs(model: cobra.Model, BBB_constituent_mets: dict) -> dict:
    """
    Get the molecular weights of the metabolites synthesized by each reaction in the passed reaction dictionary.
    Args:
        model (cobra.Model): The model to extract molecular weights from.
        BBB_constituent_mets (dict): A dictionary containing the metabolites for which to get molecular weights in the values.
    Returns:
        molecular_weights (dict): A dictionary (metabolite_ID: str, MW: float) containing the molecular weights of each metabolite.
    """
    molecular_weights = {}
    for met_ID, met_lookup in BBB_constituent_mets.items():
        MW = model.metabolites.get_by_id(met_lookup).formula_weight
        molecular_weights[met_ID] = MW
    return molecular_weights

def get_stoichiometric_coeffs(BBB_constituent_mets: dict, genmet: str, biomass_rxn: cobra.Reaction, subrxn: cobra.Reaction = None) -> dict:
    """
    Get the stoichiometric coefficients of the constituent metabolites of a biomass building block (BBB) in the biomass reaction. When there is a subreaction 
    producing a pseudo-metabolite reactant in the biomass reaction, the stoichiometric coefficients of the constituent metabolites are normalized to that of 
    the pseudo-metabolite in the biomass reaction.
    Args:
        BBB_constituent_mets (dict): A dictionary whose keys contain the metabolite IDs for which to get stoichiometric coefficients.
        genmet (str): The pseudo-metabolite ID corresponding to the BBB in the biomass reaction.
        biomass_rxn (cobra.Reaction): The biomass reaction.
        subrxn (cobra.Reaction): The subreaction for producing the pseudo-metabolite from which to take stoichiometric coefficients. Defaults to None if the 
            BBB has no associated pseudo-metabolite.        
    Returns:
        stoichiometric_coeffs (dict): A dictionary (metabolite_ID: str, stoichiometric coefficient: float) containing the stoichiometric coefficients of each 
            metabolite for the BBB.
    """
    stoichiometric_coeffs = {}
    for met_ID in BBB_constituent_mets.keys():
        if subrxn is None:
            # If there is no subreaction, simply take the stoichiometric coefficient of the metabolite in the biomass reaction.
            stoichiometric_coeffs[met_ID] = biomass_rxn.get_coefficient(met_ID)
        else:
            assert genmet in biomass_rxn.reactants, f"Pseudo-metabolite {genmet} not found in biomass reaction {biomass_rxn.id}."
            assert genmet in biomass_rxn.products, f"Pseudo-metabolite {genmet} not found in biomass synthesis subreaction {subrxn.id}."

            # take stoich coeff of constituent metabolite instead, normalized to that of the pseudo-metabolite in the biomass reaction in case it is not unity.
            stoichiometric_coeffs[met_ID] = subrxn.get_coefficient(met_ID) * (biomass_rxn.get_coefficient(genmet) / subrxn.get_coefficient(genmet))
    return stoichiometric_coeffs

def compute_BBB_mass_ratios(model: cobra.Model, biomass_rxn_id: str, BBBs_params: dict, verbose: bool = False) -> dict:
    """
    Compute the mass ratios of each biomass building block (BBB) from the biomass reaction of the model.
    Args:
        model (cobra.Model): The model to extract mass ratios from.
        biomass_rxn_id (str): The ID of the biomass reaction.
        BBBs_params (dict): A dictionary containing configuration parameters for the calculation of BBB mass ratios.
        verbose (bool): If True, prints the biomass reaction and the mass ratios. Default is False.
    Returns:
        mass_ratios (dict): A dictionary (metabolite_ID: str, mass ratio: float) containing the mass ratios of each building block.
    """
    
    assert model.reactions.has_id(biomass_rxn_id), f"Biomass reaction {biomass_rxn_id} not found in the model."    

    # Get the biomass reaction, and any subreactions (as a dict)
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)
    subreactions = {bbb:model.reactions.get_by_id(sub_rxn_id) for bbb, sub_rxn_id in BBBs_params.subrxns.items()}

    # Get the list of BBBs which are manually defined already
    manual_BBBs = [bbb for bbb in BBBs_params.BBBs if bbb not in BBBs_params.manual_mass_ratios.keys()]

    # Create a nested dictionary molecular_weights (k: BBB, v: dictionary (k: metabolite_id, v: MW))
    # molecular_weights contains the molecular weights of each constituent metabolite of a given BBB in dicts, for each BBB for which the mass ratio is NOT manually defined
    molecular_weights = {}
    # Create a similar nested dictionary for the stoichiometric coefficients in either the biomass reaction or the subreactions
    stoichiometric_coeffs = {}
    # Extract the mass ratios as a dictionary
    mass_ratios = {}

    for bbb in BBBs_params.BBBs not in manual_BBBs: # TODO: fix this line (TypeError: 'bool' object is not iterable)
        molecular_weights[bbb] = get_MWs(model, BBBs_params.constituent_mets[bbb])
        stoichiometric_coeffs[bbb] = get_stoichiometric_coeffs(BBBs_params.constituent_mets[bbb], BBBs_params.genmet[bbb], biomass_rxn, subreactions[bbb])
        # Calculate the mass ratio for each metabolite in the biomass reaction
        mass_ratios[bbb] = sum([-stoich_coeff * molecular_weights[bbb] for met, stoich_coeff in stoichiometric_coeffs[bbb].items()])/1000

    # If any mass ratios are set manually, use the values from manual_mass_ratios instead
    for bbb, mass_ratio in BBBs_params.manual_mass_ratios.items():
        assert bbb in BBBs_params.BBBs, f"Building block {bbb} not found in the BBBs in the config file."
        mass_ratios[bbb] = mass_ratio
        if verbose: print(f"Manually setting mass ratio of {bbb} to {mass_ratio} g/gDW.")

    if verbose: 
        print("Mass ratios:")
        pprint.pprint(mass_ratios)
        print(f"Sum of mass fractions: {sum(mass_ratios.values())} g/gDW.")

    return mass_ratios