import numpy as np
from scipy.optimize import fsolve

def dg_to_logk(dg_values, total_metal, total_ligand, volume_nm3, temp=300, model='sequential'):
    """
    Converts ΔG values (kJ/mol) to logK by solving population distributions.
    
    Args:
        dg_values (list): ΔG values for each step (e.g., [ΔG01, ΔG12, ...]).
        total_metal (float): Total moles of metal ions.
        total_ligand (float): Total moles of ligand.
        volume_nm3 (float): System volume in nm³.
        temp (float): Temperature in Kelvin.
        model (str): Population model type ('sequential', 'nme', etc.).
    
    Returns:
        list: Computed logK values.
    """
    R = 8.314e-3  # kJ/(mol·K)
    RT = R * temp
    n_steps = len(dg_values)
    initial_guess = [1e-5] * (n_steps + 1) + [total_ligand / (0.6022 * volume_nm3)]  # Guess for n0, n1,..., free_ligand
    
    # Define system based on model
    def equations(vars):
        populations = vars[:n_steps + 1]
        free_ligand = vars[-1]
        
        # Metal conservation
        eq_metal = total_metal / (0.6022 * volume_nm3) - sum(populations)
        
        # Ligand conservation
        ligand_used = sum((i+1)*pop for i, pop in enumerate(populations[1:]))
        eq_ligand = total_ligand / (0.6022 * volume_nm3) - free_ligand - ligand_used
        
        # Boltzmann relationships
        eqs = []
        for i in range(n_steps):
            if model == 'sequential':
                # n_{i+1} = n_i * exp(-ΔG_i / RT) * [L]
                eq = populations[i+1] - populations[i] * np.exp(-dg_values[i] / RT) * free_ligand
            elif model == 'nme':
                # Custom model (e.g., n4 depends on n3)
                # ... (extend for other models)
                pass
            eqs.append(eq)
        
        return [eq_metal, eq_ligand] + eqs
    
    solution = fsolve(equations, initial_guess)
    populations = solution[:n_steps + 1]
    free_ligand = solution[-1]
    
    # Compute logK values
    logk = []
    for i in range(n_steps):
        if populations[i] == 0:
            logk.append(np.nan)
        else:
            k = (populations[i+1] / populations[i]) / free_ligand
            logk.append(np.log10(k))
    return logk

def logk_to_dg(logk_values, total_metal, total_ligand, volume_nm3, temp=300, model='sequential'):
    """
    Converts logK values to ΔG (kJ/mol) by solving equilibrium populations.
    
    Args:
        logk_values (list): logK values for each step.
        total_metal (float): Total moles of metal ions.
        total_ligand (float): Total moles of ligand.
        volume_nm3 (float): System volume in nm³.
        temp (float): Temperature in Kelvin.
        model (str): Population model type.
    
    Returns:
        list: ΔG values (kJ/mol).
    """
    from equilibrium.conc_solver import solve_equilibrium  # Reuse equilibrium solver
    
    # Solve equilibrium to get populations and free ligand
    k_values = [10**k for k in logk_values]
    result = solve_equilibrium(total_metal, total_ligand, volume_nm3, k_values)
    populations = [result['free_metal']] + result['complexes']
    free_ligand = result['free_ligand']
    
    # Compute ΔG from populations and ligand
    R = 8.314e-3
    RT = R * temp
    dg = []
    for i in range(len(logk_values)):
        if populations[i] == 0 or free_ligand == 0:
            dg.append(np.nan)
        else:
            ratio = (populations[i+1] / (populations[i] * free_ligand))
            dg.append(-RT * np.log(ratio))
    return dg