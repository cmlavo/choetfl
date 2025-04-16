"""
Function definitions for model manipulation operations related to model setup.
"""

import cobra
import pprint
import model_io

def sanitize_model(model: cobra.Model, verbose: bool = False) -> cobra.Model:
    """
    Change all reaction and metabolite IDs that start with a digit to start with an underscore instead.
    Args:
        model (cobra.Model): The model to sanitize.
        verbose (bool): If True, prints the validation status of the model. Default is False.
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

    if verbose: print("\nModel IDs sanitized.") 

    return model

def get_MWs(model: cobra.Model, BBB_constituent_mets: dict, BBB_mw_correction: float = 0) -> dict:
    """
    Get the molecular weights of the metabolites synthesized by each reaction in the passed reaction dictionary. Adds a manually-set correction factor to the
    molecular weights of the metabolites in the model.
    Args:
        model (cobra.Model): The model to extract molecular weights from.
        BBB_constituent_mets (dict): A dictionary containing the metabolites for which to get molecular weights in the values, or 
            containing predefined molecular weights in the values to set permanently in the model.
        BBB_mw_correction (float): A dictionary containing the correction factors to add to the molecular weights of the metabolites in the model.
            Defaults to 0 if no correction is needed.
    Returns:
        molecular_weights (dict): A dictionary (metabolite_ID: str, MW: float) containing the molecular weights of each metabolite.
    """
    molecular_weights = {}
    for met_ID, met_lookup in BBB_constituent_mets.items():
        if isinstance(met_lookup, float):
            # If the value is a manually-set float, return that value
            molecular_weights[met_ID] = met_lookup + BBB_mw_correction
        elif isinstance(met_lookup, str):
            # If the value is a string, get the molecular weight from the model
            mw = model.metabolites.get_by_id(met_lookup).formula_weight + BBB_mw_correction
            assert mw is not None, f"Molecular weight of metabolite {met_lookup} not found in the model. Consider setting it manually in the config file, or changing the reference metabolite ID."
            molecular_weights[met_ID] = mw
    return molecular_weights

def get_stoichiometric_coeffs(model: cobra.Model, BBB_constituent_mets: dict, biomass_rxn: cobra.Reaction, genmet: str = None, subrxn: cobra.Reaction = None) -> dict:
    """
    Get the stoichiometric coefficients of the constituent metabolites of a biomass building block (BBB) in the biomass reaction. When there is a subreaction 
    producing a pseudo-metabolite reactant in the biomass reaction, the stoichiometric coefficients of the constituent metabolites are normalized to that of 
    the pseudo-metabolite in the biomass reaction.
    Args:
        model (cobra.Model): The model to extract stoichiometric coefficients from.
        BBB_constituent_mets (dict): A dictionary whose keys contain the metabolite IDs for which to get stoichiometric coefficients.        
        biomass_rxn (cobra.Reaction): The biomass reaction.
        genmet (str): String corresponding to the pseudo-metabolite ID of the BBB as listed in the genmets entry of the config file. Defaults to None if there 
            is no genmet.
        subrxn (cobra.Reaction): The subreaction for producing the pseudo-metabolite from which to take stoichiometric coefficients. Defaults to None if the 
            BBB has no associated pseudo-metabolite.        
    Returns:
        stoichiometric_coeffs (dict): A dictionary (metabolite_ID: str, stoichiometric coefficient: float) containing the stoichiometric coefficients of each 
            metabolite for the BBB. Note that the sign of the stiochiometric coefficient is relative to the main biomass reaction.
    """
    stoichiometric_coeffs = {}

    # Get metabolite object for genmet

    for met_ID in BBB_constituent_mets.keys():
        if subrxn is None:
            assert genmet is None, f"Pseudo-metabolite {genmet} listed in config but not associated to a subreaction."

            # If there is no subreaction, simply take the stoichiometric coefficient of the metabolite in the biomass reaction.
            stoichiometric_coeffs[met_ID] = biomass_rxn.get_coefficient(met_ID)
        else:
            genmet_obj = model.metabolites.get_by_id(genmet)

            assert genmet_obj in biomass_rxn.reactants, f"Pseudo-metabolite {genmet_obj.id} not found in biomass reaction {biomass_rxn.id}."
            assert genmet_obj in subrxn.products, f"Pseudo-metabolite {genmet_obj.id} not found in biomass synthesis subreaction {subrxn.id}."

            # take stoich coeff of constituent metabolite instead, normalized to that of the pseudo-metabolite in the biomass reaction in case it is not unity.
            stoichiometric_coeffs[met_ID] = subrxn.get_coefficient(met_ID) * -(biomass_rxn.get_coefficient(genmet) / subrxn.get_coefficient(genmet))
     
    return stoichiometric_coeffs

def compute_BBB_mass_ratios(model: cobra.Model, biomass_rxn_id: str, BBBs_params: dict, verbose: bool = False) -> dict:
    """
    Compute the mass ratios of each biomass building block (BBB) from the biomass reaction of the model. Mass ratios set manually in the
    configuration file overwrite the calculated mass ratios.
    Args:
        model (cobra.Model): The model to extract mass ratios from.
        biomass_rxn_id (str): The ID of the biomass reaction.
        BBBs_params (dict): A dictionary containing configuration parameters for the calculation of BBB mass ratios.
        verbose (bool): If True, prints mass ratios and status updates. Default is False.
    Returns:
        mass_ratios (dict): A dictionary (metabolite_ID: str, mass ratio: float) containing the mass ratios of each building block.
    """
    
    assert model.reactions.has_id(biomass_rxn_id), f"Biomass reaction {biomass_rxn_id} not found in the model."    

    # Get the biomass reaction, and any subreactions (as a dict)
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)
    subreactions = {bbb:model.reactions.get_by_id(sub_rxn_id) for bbb, sub_rxn_id in BBBs_params.subrxns.items()}

    # Extract the mass ratios as a dictionary
    mass_ratios = {}

    if verbose: print("\n")
    
    for bbb in BBBs_params.BBBs:

        if bbb in BBBs_params.manual_mass_ratios.keys():
            # mass ratio already manually defined
            mass_ratios[bbb] = BBBs_params.manual_mass_ratios.get(bbb)
            if verbose: print(f"Manually set mass ratio of {bbb} to {mass_ratios[bbb]:.5f} g/gDW.")

        else:
            # compute mass ratio
            molecular_weights = get_MWs(model, BBBs_params.constituent_mets[bbb], BBBs_params.mw_correction.get(bbb, 0))
            stoichiometric_coeffs = get_stoichiometric_coeffs(model, BBBs_params.constituent_mets[bbb], biomass_rxn, BBBs_params.genmets.get(bbb, None), subreactions.get(bbb, None))
            #breakpoint()
            # Calculate the mass ratio for each metabolite in the biomass reaction
            mass_ratios[bbb] = sum([-stoich_coeff * molecular_weights[met] for met, stoich_coeff in stoichiometric_coeffs.items()])/1000
            if verbose: print(f"Computed mass ratio of {bbb} as {mass_ratios[bbb]:.5f} g/gDW.")

    if verbose: 
        print("Mass ratios:")
        pprint.pprint(mass_ratios)
        print(f"Sum of mass fractions: {sum(mass_ratios.values())} g/gDW.")

    if BBBs_params.mass_ratio_normalization:
        # Normalize the mass ratios to sum to 1
        total_mass = sum(mass_ratios.values())
        for bbb in mass_ratios.keys():
            mass_ratios[bbb] /= total_mass
        if verbose:
            print("Normalized mass ratios:")
            pprint.pprint(mass_ratios)

    return mass_ratios

def modify_GAM(model: cobra.Model, biomass_rxn_id: str, BBBs_params: dict, GAM_params: dict, verbose: bool = False) -> cobra.Model:
    """
    Modify the growth-associated maintenance (GAM) reaction of the model to remove the energetic cost of growth, maintenance
    and peptide polymerization from the model. The energetic cost will later be added explicitly, so removing it from
    the GAM avoids overestimating this energetic requirement.
    Args:
        model (cobra.Model): The model to modify.
        biomass_rxn_id (str): The ID of the biomass reaction.
        BBBs_params (dict): A dictionary containing configuration parameters for the BBBs.
        GAM_params (dict): A dictionary containing configuration parameters for the GAM modification.
        verbose (bool): If True, prints status updates. Default is False.
    Returns:
        GAM_modified_model (cobra.Model): The modified model.
    """
    
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)

    subrxn_id = BBBs_params.subrxns.get('protein', None)
    subrxn = None if subrxn_id == None else model.reactions.get_by_id(subrxn_id)

    # getting number of moles of amino acids required to synthesize 1gDW of biomass
    aa_stoichiometric_coeffs = get_stoichiometric_coeffs(model, BBBs_params.constituent_mets['protein'], biomass_rxn, BBBs_params.genmets.get('protein', None), subrxn)
    total_moles_aa = -sum(aa_stoichiometric_coeffs.values()) # negative value since they are in the reactants, *-1  to make it positive

    # since 2 moles of GTP are required to attach 1 mole of aa to the growing peptide chain, the total GTP expense is 2 * total_moles_aa
    gtp_expense = 2 * total_moles_aa

    # create dict of GAM metabolite (cobra.Metabolite) : GTP expense (float) pairs to subtract these values from the biomass reaction
    GAM_metabolites = [model.metabolites.get_by_id(met_ID) for met_ID in GAM_params.GAM_mets]
    GAM_dict = {}
    for met in GAM_metabolites: 
        assert met in biomass_rxn.metabolites, f"GAM metabolite {met.id} not found in biomass reaction {biomass_rxn.id}. Modify etfl.GAM_params.GAM_mets entry in config file."
        GAM_dict[met] = -gtp_expense if met in biomass_rxn.reactants else gtp_expense

    # subtract the GTP expense from the biomass reaction
    biomass_rxn.subtract_metabolites(GAM_dict)

    if verbose: 
        print("\n")
        print(f"GTP expense for peptide polymerization: {gtp_expense:.5f} mmol/gDW")
        print("Modified biomass reaction:")
        print(biomass_rxn.reaction)

    return model

def develop_coupling_dict(model: cobra.Model, enz_coupl_params: dict, verbose: bool = False) -> dict:
    """
    Develop a dictionary coupling all metabolic reactions to a list of etfl.Enzyme objects.
    Args:
        model (cobra.Model): The model to extract coupling information from.
        enz_coupl_params (dict): Configuration parameters for the creation of the coupling dictionary.
        verbose (bool): If True, prints status updates. Default is False.
    Returns:
        coupling_dict (dict): A dictionary (cobra.Reaction: cobra.Enzyme) containing the coupling information for each reaction in the model.
    """
    
    # Importing protein complex data as dataframes from the .csv files
    prot_complx_dfs = [model_io.import_protein_complex_data(file_name, verbose) for file_name in enz_coupl_params.prot_complx_files]

    # Merging the dataframes into a single dataframe
    prot_complx_df = model_io.merge_protein_complex_data(prot_complx_dfs, enz_coupl_params.comp_coeff_mismatch_treatment, verbose)
    breakpoint()
