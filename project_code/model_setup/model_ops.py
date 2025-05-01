"""
Function definitions for model manipulation operations related to model setup.
"""

import cobra
from ordered_set import OrderedSet
import pandas as pd
from pprint import pprint
import model_io
import utils
from etfl.core import ThermoMEModel

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

    if verbose: print("")
    
    for bbb in BBBs_params.BBBs:

        if bbb in BBBs_params.manual_mass_ratios.keys():
            # mass ratio already manually defined
            mass_ratios[bbb] = BBBs_params.manual_mass_ratios.get(bbb)
            if verbose: print(f"Manually set mass ratio of {bbb} to {mass_ratios[bbb]:.5f} g/gDW.")

        else:
            # compute mass ratio
            molecular_weights = get_MWs(model, BBBs_params.constituent_mets[bbb], BBBs_params.mw_correction.get(bbb, 0))
            stoichiometric_coeffs = get_stoichiometric_coeffs(model, BBBs_params.constituent_mets[bbb], biomass_rxn, BBBs_params.genmets.get(bbb, None), subreactions.get(bbb, None))
            # Calculate the mass ratio for each metabolite in the biomass reaction
            mass_ratios[bbb] = sum([-stoich_coeff * molecular_weights[met] for met, stoich_coeff in stoichiometric_coeffs.items()])/1000
            if verbose: print(f"Computed mass ratio of {bbb} as {mass_ratios[bbb]:.5f} g/gDW.")

    if verbose: 
        print("Mass ratios:")
        pprint(mass_ratios)
        print(f"Sum of mass fractions: {sum(mass_ratios.values())} g/gDW.")

    if BBBs_params.mass_ratio_normalization:
        # Normalize the mass ratios to sum to 1
        total_mass = sum(mass_ratios.values())
        for bbb in mass_ratios.keys():
            mass_ratios[bbb] /= total_mass
        if verbose:
            print("Normalized mass ratios:")
            pprint(mass_ratios)

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
        print(f"\nGTP expense for peptide polymerization: {gtp_expense:.5f} mmol/gDW")
        print("Modified biomass reaction:")
        print(biomass_rxn.reaction)

    return model

def assoc_rxn_enz(model: cobra.Model, enzyme_dict: dict, genes_dict: dict, enz_coupl_params: dict, verbose: bool = False) -> dict:
    """
    Associate reactions in the model with their corresponding enzyme(s) based on the provided model, enzyme dictionary, and gene dictionary.
    Set the appropriate k values for the enzymes based on the kcat data and the configuration file.
    Args:
        model (cobra.Model): The cobra model.
        enzyme_dict (dict): A dictionary containing enzyme IDs as keys and etfl.Enzyme objects as values.
        genes_dict (dict): A dictionary containing enzyme IDs as keys and OrderedSet of gene IDs as values.
        enz_coupl_params (dict): Configuration parameters for the enzyme coupling.
        verbose (bool): If True, prints status updates. Default is False.
    Returns:
        rxn_enz_dict (dict): A dictionary containing cobra.Reaction objects as keys and lists of etfl.Enzyme objects, with correct k values, as values.
    """
    
    rxn_enz_dict = {}
    init_num_enzymes = len(enzyme_dict)

    kcat_df = model_io.import_kcat_data(enz_coupl_params.kcat_file, enz_coupl_params.k_cat_default, verbose)

    # Loop through each reaction in the model
    for rxn in model.reactions:
        # Get the list of genes related to the reaction in the model
        rxn_genes = OrderedSet([gene.id for gene in rxn.genes])

        # Compute the Jaccard similarity between the genes for the reaction, and the genes for each enzyme, to try to associate the reaction with an enzyme. Store in a dataframe.
        jaccard_scores = pd.Series({enz: utils.jaccard_score(rxn_genes, enz_genes) for enz, enz_genes in genes_dict.items()})
        matched_enz = [enzyme_dict[score.index] for score in jaccard_scores if abs(jaccard_scores - 1.) <= enz_coupl_params.jaccard_tol] # List of etfl.Enzyme objects which are a perfect Jaccard match

        # If there are no perfect matches
        # TODO: add further logic to catch non-perfect matches and add them to matched_enz. For example:
        #   specific reactions which need manual enzyme assignment
        #   for reactions with two enzymes but only one is known, infer composition of the other
        #   for genes missing from the GEM model but present in the protein complex data
        #   monomeric enzymes (only one gene in the gene set), removing those which are part of a larger complex to avoid redunancy
        #   infer from GPR
        # Only if enzyme_inference is True, else do not infer any data and only use what is given.
        ...
        
        if len(matched_enz) > 0:
            rxn_enz_dict[rxn] = list(set(matched_enz)) # ensure that the list of enzymes is unique

        else:
            # TODO: add logic for various cases where no enzyme(s) are found for a given reaction. For example:
            #   Add monomeric enzymes for each gene assigned to an unmatched reaction.
            #   For enzymes with multiple kcat values, create new enzyme for each individual kcat
            # Add new enzymes to enzyme_dict (only add new enzymes if enzyme_inference is True)
            ...

        for enz in rxn_enz_dict[rxn]:
            # find enzyme in the dataframe
            row = kcat_df[kcat_df['prot_id'] == enz.id]

            # If the enzyme is found in the dataframe, extract the kcat value(s) from the dataframe and set it in the Enzyme object, else leave it as the default value
            if not row.empty: 
                enz.kcat_fwd = row['k_cat_fwd']
                enz.kcat_bwd = row['k_cat_bwd'] 

        # TODO: add kcat fwd/bwd values depending on the logic (monomeric enzymes, transport reactions, manual values etc.)
        ...
        
    if verbose:
        if enz_coupl_params.enzyme_inference: print(f"Inference step created {len(enzyme_dict) - init_num_enzymes} new enzymes.")
        print(f"{len(enzyme_dict)} enzymes successfully associated to {len(rxn_enz_dict)} reactions.")

    return rxn_enz_dict

def develop_coupling_dict(model: cobra.Model, enz_coupl_params: dict, verbose: bool = False) -> dict:
    """
    Develop a dictionary coupling several metabolic reactions to a list of etfl.Enzyme objects.
    Args:
        model (cobra.Model): The model to extract coupling information from.
        enz_coupl_params (dict): Configuration parameters for the creation of the coupling dictionary.
        verbose (bool): If True, prints status updates. Default is False.
    Returns:
        coupling_dict (dict): A dictionary (cobra.Reaction: list of cobra.Enzyme) containing the coupling information for some reactions in the model.
    """

    if verbose: print("")
    
    # Importing protein complex data as dataframes from the .csv files
    prot_complx_dfs = [model_io.import_protein_complex_data(file_name, verbose) for file_name in enz_coupl_params.prot_complx_files]

    if verbose: print("")

    # Merging the dataframes into a single dataframe
    prot_complx_df = model_io.merge_protein_complex_data(prot_complx_dfs, enz_coupl_params.comp_coeff_mismatch_treatment, verbose)

    if verbose: print("")

    # Creating a dictionary of protein-complex:etfl.Enzyme objects from the protein complex dataframe
    enzyme_dict = model_io.assoc_prot_enz(prot_complx_df, verbose)

    # Creating a dictionary of protein-complex (str) : genes ID (OrderedSet) pairs from the protein complex dataframe
    genes_dict = {prot_ID : OrderedSet(genes) for prot_ID, genes in zip(prot_complx_df['prot_id'], prot_complx_df['gene_id'])}

    if verbose: print("")

    # Couple as many cobra reactions as possible to to one or multiple catalyzing enzymes (in a dict)
    rxn_enz_dict = assoc_rxn_enz(model, enzyme_dict, genes_dict, enz_coupl_params, verbose)
    
    return rxn_enz_dict

def create_ETFL_model(model: cobra.Model, baseline_config: dict, etfl_config: dict, solver_config: dict, verbose: bool = False) -> ThermoMEModel:
    """
    Create the ETFL model from the baseline model, loading in thermo data from the data folder specified in the config file
    and setting solver configuration parameters.
    Args:
        model (cobra.Model): The baseline model to use as a template for the ETFL model.
        baseline_config (dict): Configuration parameters for the baseline model.
        etfl_config (dict): Configuration parameters for the ETFL model.
        solver_config (dict): Configuration parameters for the solver.
        verbose (bool): If True, prints status updates. Default is False.
    Returns:
        ETFL_model (ThermoMEModel): The created ETFL model.
    """

    if verbose: print("")

    # Load thermodynamic database
    thermo_data = model_io.import_thermo_data(etfl_config.thermo_data.file_name, verbose)

    if verbose: print("")    

    # Create the ETFL model
    ETFL_model = ThermoMEModel(thermo_data = thermo_data, 
                               model = model, 
                               name = etfl_config.model_name,
                               growth_reaction = baseline_config.biomass_rxn_id,
                               mu_range = (etfl_config.mu_discr.mu_min, etfl_config.mu_discr.mu_max),
                               n_mu_bins = etfl_config.mu_discr.mu_bins,
                               solver = solver_config.name,                               
                               )
    
    # Additional model setup
    ETFL_model.solver.configuration.verbosity = int(verbose)
    ETFL_model.solver.configuration.tolerances.feasibility = solver_config.tolerance
    ETFL_model.solver.configuration.timeout = solver_config.timeout

    # Prepare and convert the model for TFBA
    ETFL_model.prepare()
    ETFL_model.convert()

    if verbose: print(f"Sucessfully created ETFL model {etfl_config.name}.")

    return ETFL_model

def fix_protein_alloc(model: ThermoMEModel, protein_mass_ratio: float) -> None:
    ...

def fix_RNA_alloc(model: ThermoMEModel, rna_mass_ratio: float) -> None:
    ...

def fix_DNA_alloc(model: ThermoMEModel, dna_mass_ratio: float, gc_ratio: float, chromosome_length: float) -> None:
    ...

def constrain_enzymes(model: ThermoMEModel, protein_mass_ratio: dict) -> None:
    ...

def update_copy_numbers(model: ThermoMEModel, copy_dict: dict) -> None:
    ...

def develop_trna_enz_coupling_dict(model: ThermoMEModel, charging_enz_dict: dict) -> dict:
    ...

def populate_ETFL_model(model: ThermoMEModel, mass_ratios: dict, rxn_enz_coupling_dict: dict, etfl_config: dict) -> None:
    """
    Populate the ETFL model with proteomics, metabolomics, and transcriptomics data. Applies modifications directly to the passed
    model, so has no return value.
    Args:
        model (ThermoMEModel): The ETFL model to populate.
        mass_ratios (dict): A dictionary containing the mass ratios of the biomass building blocks.
        rxn_enz_coupling_dict (dict): A dictionary containing the coupling information for the reactions in the model.
        etfl_config (dict): Configuration parameters for the ETFL model.
    """

    nucleotides_sequences = model_io.get_nucleotides_sequences(...)
    model.add_nucleotide_sequences(nucleotides_sequences)

    aa_sequences = model_io.get_amino_acid_sequences(...)
    model.add_peptide_sequences(aa_sequences)
    
    essential_BBs = model_io.get_essential_BBs(...)
    aa_dict = model_io.get_aa_dict(...)
    rna_nucleotides_tp = model_io.get_rna_nucleotides_tp(...)
    rna_nucleotides_mp = model_io.get_rna_nucleotides_mp(...)
    model.add_essentials(essential_BBs, aa_dict, rna_nucleotides_tp, rna_nucleotides_mp)

    mrna_dict = model_io.get_mrna_dict(...)
    model.add_mrnas(mrna_dict)

    for rib_ID, free_rib_ratio in etfl_config.ribosome_ratios:
        ribosome = model_io.get_ribosomes(rib_ID)
        model.add_ribosome(ribosome, free_rib_ratio)

    for rnap_ID, free_rnap_ratio in etfl_config.rnap_ratios:
        rnap = model_io.get_rnap(rnap_ID)
        model.add_rnap(rnap, free_rnap_ratio)

    transcription_dict = model_io.get_transcription_dict(...)
    model.add_transcription_by(transcription_dict)

    translation_dict = model_io.get_translation_dict(...)
    model.add_translation_by(translation_dict)

    model.build_expression()
    model.add_enzymatic_coupling(rxn_enz_coupling_dict)

    charging_enz_list, charging_enz_dict = model_io.get_trna_charging_enzymes(...)
    model.add_enzymes(charging_enz_list)

    model.add_dummies(nt_ratios = etfl_config.nt_ratios, 
                      mrna_kdeg = etfl_config.k_constants.k_deg_mrna, 
                      mrna_length = etfl_config.mrna_length_avg, 
                      aa_ratios = etfl_config.aa_ratios, 
                      enzyme_kdeg = etfl_config.k_constants.k_deg_enz, 
                      peptide_length = etfl_config.peptide_length_avg)
    
    if etfl_config.variable_allocation:
        # TODO: add logic for variable allocation
        ...
    else:
        fix_protein_alloc(model, mass_ratios['protein'])
        fix_RNA_alloc(model, mass_ratios['RNA'])
        fix_DNA_alloc(model, mass_ratios['DNA'], etfl_config.gc_ratio, etfl_config.chromosome_length)

    update_copy_numbers(model, etfl_config.copy_dict)

    model.populate_expression()

    model.add_trna_mass_balances()

    trna_enz_coupling_dict = develop_trna_enz_coupling_dict(model, charging_enz_dict)
    model.add_enzymatic_coupling(trna_enz_coupling_dict)

    if etfl_config.constrain_enzymes: constrain_enzymes(model, mass_ratios['protein'])

    model.repair()