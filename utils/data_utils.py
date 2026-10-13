import pandas as pd
import numpy as np
import yaml
from yaml.nodes import ScalarNode, SequenceNode, MappingNode
import cobra
from cobra.flux_analysis import pfba
import warnings
import matplotlib.pyplot as plt
from tqdm import tqdm


#Define custom loader to handle MATLAB's empty !!omap scalars
class MatlabYamlSafeLoader(yaml.SafeLoader):
    pass

def custom_omap_constructor(loader, node):
    """Safely handle MATLAB/RAVEN's malformed empty !!omap scalars."""
    if isinstance(node, ScalarNode):
        return {}  # Return empty dictionary if it's an empty scalar
    elif isinstance(node, MappingNode):
        return loader.construct_mapping(node)
    elif isinstance(node, SequenceNode):
        # Normal omap: list of single-key dicts -> combine into a single dict
        pairs = []
        for subnode in node.value:
            pairs.extend(loader.construct_pairs(subnode))
        return dict(pairs)
    return {}

# Register the constructor for standard YAML omap tags
MatlabYamlSafeLoader.add_constructor('tag:yaml.org,2002:omap', custom_omap_constructor)
MatlabYamlSafeLoader.add_constructor('!omap', custom_omap_constructor)


#Categorise boundary reactions into imports and exports, accounting for stoichiometry and pseudo-metabolites
def analyze_boundary_reactions(model, show_pseudometabolites=False):
    """
    Robustly identifies and categorizes all boundary/exchange reactions into
    Imports (uptakes) and Exports (effluxes).
    """
    records = []
    
    for rxn in model.boundary:
        # Boundary reactions have exactly 1 metabolite
        met = list(rxn.metabolites.keys())[0]
        coeff = rxn.metabolites[met]
        
        # Check if it's an enzyme/protein pseudo-reaction (GECKO specific)
        is_pseudo = any(term in rxn.id.lower() or term in met.id.lower() 
                        for term in ["prot_", "pool", "draw_prot"])
        
        if not show_pseudometabolites and is_pseudo:
            continue
            
        # Determine actual import/export capacity accounting for stoichiometry
        # Standard convention: met[e] <=>  (coeff == -1)
        if coeff < 0:
            can_import = rxn.lower_bound < 0
            can_export = rxn.upper_bound > 0
            max_import = -rxn.lower_bound if can_import else 0.0
            max_export = rxn.upper_bound if can_export else 0.0
        # Inverted convention: <=> met[e] (coeff == +1)
        else:
            can_import = rxn.upper_bound > 0
            can_export = rxn.lower_bound < 0
            max_import = rxn.upper_bound if can_import else 0.0
            max_export = -rxn.lower_bound if can_export else 0.0
            
        records.append({
            "Reaction_ID": rxn.id,
            "Reaction_Name": rxn.name,
            "Metabolite_ID": met.id,
            "Metabolite_Name": met.name,
            "Compartment": met.compartment,
            "Lower_Bound": rxn.lower_bound,
            "Upper_Bound": rxn.upper_bound,
            "Allows_Import": can_import,
            "Max_Import_Rate": max_import,
            "Allows_Export": can_export,
            "Max_Export_Rate": max_export,
            "Is_Pseudo": is_pseudo
        })
        
    df = pd.DataFrame(records)
    return df

#proteome allocation analysis for GECKO models
def get_proteome_allocation(model, target_rxn="1687"):
    m = model.copy()
    m.objective = target_rxn
    sol = m.optimize()
    cap = abs(m.reactions.get_by_id("prot_pool_exchange").lower_bound)

    records = []
    for rxn in m.reactions:
        # 1. Skip the pool delivery/exchange helper reactions
        if rxn.id.startswith("usage_prot_") or rxn.id == "prot_pool_exchange":
            continue
            
        for met, coeff in rxn.metabolites.items():
            # 2. Only match specific enzymes, NOT the global prot_pool metabolite
            if met.id.startswith("prot_") and met.id != "prot_pool" and coeff < 0:
                mass = abs(sol.fluxes.get(rxn.id, 0.0) * coeff)
                if mass > 1e-6:
                    sub = rxn.subsystem
                    sub_str = ", ".join(sub) if isinstance(sub, list) else (sub or "Other")
                    records.append({"subsystem": sub_str, "mass_mg": mass})

    df = pd.DataFrame(records).groupby("subsystem", as_index=False)["mass_mg"].sum()
    df["%_pool"] = (df["mass_mg"] / cap) * 100
    return df.sort_values(by="mass_mg", ascending=False).reset_index(drop=True)

#FSEOF implementation based on De Biaggi et al. (2026) workflow, with fixes for MATLAB quirks and solver differences.
def run_fseof_debiaggi(
    model, 
    prod_target_rxn, 
    bio_rxn="xBIOMASS", 
    n_steps=16,
    filter_complex_rxns=True,       # Driver passes true for every call
    use_gurobi_if_available=True
):
    """
    Exact replication of De Biaggi et al. (2026) FSEOFalt workflow.
    """
  
    m = model.copy()

    # unfortunately dont have access to the original solver configuration, so we will just use the default solver.
    print(f"Using solver: {m.solver.interface.__name__}")
    m.reactions.get_by_id(prod_target_rxn).lower_bound = 0
    m.reactions.get_by_id(prod_target_rxn).upper_bound = 1000.0

    # --- FIX 4: Gene-rule filtering (String-only check, matching MATLAB withGR) ---
    # MATLAB: withGR = ~cellfun(@isempty, model.grRules); stdIdx = contains(model.grRules,'standard'); withGR(stdIdx)=0;
    valid_rxn_ids = []
    for r in m.reactions:
        rule = r.gene_reaction_rule
        if rule and rule.strip() and ("standard" not in rule.lower()):
            valid_rxn_ids.append(r.id)

    # Initial target production at max biomass using standard FBA
    # MATLAB: model = setParam(model, 'obj', params.bioRxn, 1); sol = solveLP(model);
    m.objective = bio_rxn
    m.reactions.get_by_id(bio_rxn).lower_bound = 0.0
    sol_ini = m.optimize()
    ini_target = sol_ini.fluxes[prod_target_rxn]

    # Theoretical maximum target yield using pFBA
    # MATLAB: model = setParam(model,'obj',prodTargetRxn,1); sol = solveLP(model,1);
    m.objective = prod_target_rxn
    try:
        sol_max = pfba(m)
    except Exception:
        sol_max = m.optimize()
    max_target = sol_max.fluxes[prod_target_rxn] * 0.9  # 90% threshold

    print(f"Scanning target '{prod_target_rxn}' from {ini_target:.6e} to {max_target:.6e} ({n_steps} steps)")

    # Define alpha grid and flux matrix
    alpha = np.linspace(ini_target, max_target, num=n_steps)
    v_matrix = np.zeros((len(m.reactions), n_steps))

    # Iterative FSEOF Scan
    for i, val in enumerate(alpha):
        with m:
            # MATLAB: model = setParam(model, 'eq', prodTargetRxnIdx, alpha(i));
            m.reactions.get_by_id(prod_target_rxn).bounds = (val, val)
            
            # MATLAB: model.lb(bioRxnIdx) = 0; model = setParam(model, 'obj', params.bioRxn, 1);
            m.reactions.get_by_id(bio_rxn).lower_bound = 0.0
            m.objective = bio_rxn
            
            # MATLAB: sol = solveLP(model, 1);
            try:
                sol_step = pfba(m)
                v_matrix[:, i] = [sol_step.fluxes[r.id] for r in m.reactions]
            except Exception:
                v_matrix[:, i] = 0.0

    all_rxn_ids = [r.id for r in m.reactions]
    df_v = pd.DataFrame(v_matrix, index=all_rxn_ids).loc[valid_rxn_ids]

    # Filter out always-zero reactions across all steps
    # MATLAB: zero_flux = ~all(abs(v_matrix(:,1:nSteps))==0,2);
    df_v = df_v.loc[(np.abs(df_v.values) > 1e-9).any(axis=1)]

    # --- FIX 2: filterComplexRxns (Applied BEFORE slope calculation & quantile) ---
    # MATLAB driver explicitly passed: filterComplexRxns = true
    # MATLAB: complex_rxns = contains(model.grRules(...), 'and'); target_rxns = target_rxns(~complex_rxns);
    if filter_complex_rxns:
        non_complex_ids = []
        for r_id in df_v.index:
            rule = m.reactions.get_by_id(r_id).gene_reaction_rule.lower()
            # Match MATLAB's contains(..., 'and')
            if ' and ' not in rule and not rule.startswith('and ') and not rule.endswith(' and'):
                non_complex_ids.append(r_id)
        
        n_dropped = len(df_v) - len(non_complex_ids)
        df_v = df_v.loc[non_complex_ids]
        print(f"Filtered out {n_dropped} complex reactions containing 'and' (filterComplexRxns=True).")

    # Slope denominator quirk from MATLAB: abs(maxTarget - maxTarget/nSteps - 1)
    matlab_denominator = abs(max_target - (max_target / n_steps) - 1.0)

    targets = []
    for rxn_id in df_v.index:
        raw_fluxes = df_v.loc[rxn_id].values
        abs_fluxes = np.abs(raw_fluxes)

        # Monotonicity checks
        # If using Gurobi, clean machine epsilon rounding is used; 
        # if using GLPK, allow for solver plateaus
        diffs = np.diff(abs_fluxes)
        
        if m.solver.interface.__name__ == 'optlang.gurobi_interface':
            is_ascend = np.all(diffs > 1e-12)
            is_descend = np.all(diffs < -1e-12)
        else:
            # GLPK fallback: allow flat steps from degeneracy, require net movement
            is_ascend = np.all(diffs >= -1e-9) and (abs_fluxes[-1] - abs_fluxes[0] > 1e-6)
            is_descend = np.all(diffs <= 1e-9) and (abs_fluxes[0] - abs_fluxes[-1] > 1e-6)

        if is_ascend or is_descend:
            # MATLAB: abs(v_matrix(i, nSteps-1) - v_matrix(i, 1)) / denominator
            slope = abs(abs_fluxes[-2] - abs_fluxes[0]) / matlab_denominator
            
            if is_ascend:
                action_type = "OE"
            else:
                action_type = "KO" if abs(raw_fluxes[-1]) < 1e-9 else "KD"

            targets.append({
                'rxn_id': rxn_id,
                'name': m.reactions.get_by_id(rxn_id).name,
                'slope': slope,
                'action': action_type,
                'grRule': m.reactions.get_by_id(rxn_id).gene_reaction_rule
            })

    df_targets = pd.DataFrame(targets)
    if df_targets.empty:
        print("No monotonic targets found.")
        return pd.DataFrame()

    # Sort descending by slope
    df_targets = df_targets.sort_values(by='slope', ascending=False).reset_index(drop=True)
    n_before_q75 = len(df_targets)

    # --- FIX 6: Quantile Calculation Matching MATLAB ---
    # MATLAB quantile() uses midpoint/type-5 interpolation by default, 
    # whereas pandas defaults to 'linear'. 'midpoint' is closer.
    q75 = df_targets['slope'].quantile(0.75, interpolation='midpoint')
    
    # MATLAB: non_zero_slope = slope_rxns > quantile(slope_rxns, 0.75);
    df_filtered = df_targets[df_targets['slope'] > q75].copy()

    # Remove any stray usage_prot_ reactions
    # MATLAB: toKeep = ~startsWith(target_rxns, 'usage_prot_');
    df_filtered = df_filtered[~df_filtered['rxn_id'].str.startswith('usage_prot_')].reset_index(drop=True)

    print(f"Identified {n_before_q75} monotonic reactions (before 75% cutoff).")
    print(f"Calculated 75th percentile cutoff: {q75:.6f}")
    print(f"Retained {len(df_filtered)} reactions after combined 75th percentile filter.")
    print(f"   -> OE Reactions : {len(df_filtered[df_filtered['action']=='OE'])}")
    print(f"   -> KD/KO Reactions: {len(df_filtered[df_filtered['action'].isin(['KD','KO'])])}")

    return df_filtered


#FCC calculation for GECKO models, based on the authors' MATLAB FCC.m script.
def calculate_fcc_kcat(model, target_rxn, Qa=1.001, top_n=10):
    """
    Direct Python translation of the authors' MATLAB FCC.m script.
    Perturbs kcat by Qa (+0.1%) for each enzyme-associated reaction and 
    computes normalized Flux Control Coefficients (FCCs).
    
    Parameters:
    - model (cobra.Model): A GECKO _pooled.yml model (e.g. ecOBEe_pooled.yml)
    - target_rxn (str): Objective reaction ID (e.g. 'EXC_OUT_m1640')
    - Qa (float): Perturbation factor (default 1.001 = +0.1%)
    - top_n (int): Number of top FCCs to return and plot
    
    Returns:
    - df_fcc (pd.DataFrame): Sorted table with FCCs, reactions, and enzymes
    """
    m = model.copy()
    m.reactions.get_by_id(target_rxn).bounds = (0.0, 1000.0)
    # 1. Set target as objective & calculate reference target flux
    # (Matches MATLAB: ecModel = setParam(ecModel, 'obj', modelFCC.obj{i}, 1); sol = solveLP(ecModel);)
    m.objective = target_rxn
    sol_ref = m.optimize()
    target_ref = sol_ref.objective_value
    
    if target_ref is None or target_ref <= 1e-9:
        raise ValueError(f"Baseline flux for {target_rxn} is 0 or infeasible. Check model constraints.")

    print(f"[{m.id if m.id else 'Model'}] Target: {target_rxn} | Baseline Flux: {target_ref:.6e} mmol/gDW/h")

    # 2. Identify all enzyme-catalyzed reactions and their protein metabolites
    # In GECKO ecModels, enzyme-catalyzed reactions consume 'prot_' pseudo-metabolites
    kcat_targets = []
    for rxn in m.reactions:
        if not rxn.id.startswith("usage_prot_") and rxn.id != "prot_pool_exchange":
            for met, coeff in rxn.metabolites.items():
                if met.id.startswith("prot_") and coeff < 0:
                    enz_id = met.id.replace("prot_", "")
                    kcat_targets.append({
                        'rxn_id': rxn.id,
                        'rxn_name': rxn.name,
                        'met': met,
                        'enz_id': enz_id,
                        'orig_coeff': coeff
                    })

    print(f"Identified {len(kcat_targets)} enzyme-associated reaction steps. Calculating FCCs...")

    # 3. Iteratively perturb kcat and calculate FCC
    # (Matches MATLAB loop lines 56-70)
    raw_fccs = []
    
    for item in tqdm(kcat_targets, desc="Calculating FCCs"):
        rxn = m.reactions.get_by_id(item['rxn_id'])
        prot_met = item['met']
        orig_coeff = item['orig_coeff']
        
        # In GECKO, coeff = -MW/kcat.
        # Multiplying kcat by Qa is mathematically identical to dividing coeff by Qa:
        new_coeff = orig_coeff / Qa
        delta_coeff = new_coeff - orig_coeff
        
        # Apply perturbation
        rxn.add_metabolites({prot_met: delta_coeff})
        
        # Solve perturbed model (Matches MATLAB: sol = solveLP(ecModel, 1);)
        sol_pert = m.optimize()
        target_fcc = sol_pert.objective_value if sol_pert.status == "optimal" else target_ref
        
        # Calculate FCC: (targetFCC - targetRef) / (targetRef * (Qa - 1))
        # (Matches MATLAB line 65)
        fcc_val = (target_fcc - target_ref) / (target_ref * (Qa - 1.0))
        raw_fccs.append(max(0.0, fcc_val))  # Ignore negative numerical noise
        
        # Restore original coefficient
        rxn.add_metabolites({prot_met: -delta_coeff})

    # 4. Build output DataFrame & Normalize
    # (Matches MATLAB line 76: result.FCCs = FCCs/sum(FCCs);)
    df_fcc = pd.DataFrame(kcat_targets)
    df_fcc['raw_FCC'] = raw_fccs
    
    total_fcc = df_fcc['raw_FCC'].sum()
    if total_fcc > 1e-12:
        df_fcc['FCC'] = df_fcc['raw_FCC'] / total_fcc
    else:
        df_fcc['FCC'] = 0.0

    # 5. Sort descending
    df_fcc = df_fcc.sort_values(by='FCC', ascending=False).reset_index(drop=True)
    
    # Drop internal metabolite object for clean display
    df_fcc = df_fcc.drop(columns=['met', 'orig_coeff'])
    
    return df_fcc