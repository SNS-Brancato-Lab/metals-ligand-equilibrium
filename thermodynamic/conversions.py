import numpy as np
from scipy.optimize import fsolve

def dg_to_logk(dg_values, total_metal, total_ligand, volume_nm3, temp=300, model='sequential', dg_std = None):
    """
    Converts ΔG values (kJ/mol) to logK by solving population distributions.
    
    Args:
        dg_values (list): ΔG values for each step (e.g., [ΔG01, ΔG12, ...]).
        total_metal (float): Total moles of metal ions.
        total_ligand (float): Total moles of ligand.
        volume_nm3 (float): System volume in nm³.
        temp (float): Temperature in Kelvin.
        model (str): Population model type ('sequential', 'nme', etc.).
        dg_sdt (list): Standard distribution of each ΔG values.  
    
    Returns:
        list: Computed logK values.

    Raises:
        TypeError: Raised when dg_std is not correctly provided.
    """
    R = 8.314e-3  # kJ/(mol·K)
    RT = R * temp
    n_steps = len(dg_values)
    initial_guess = [1e-5] * (n_steps + 1) + [total_ligand / (0.6022 * volume_nm3)]  # Guess for n0, n1,..., free_ligand
    
    # Define system based on model
    def equations(tot_populations, *args):
        
        populations = tot_populations[:n_steps + 1]
        free_ligand = tot_populations[-1]
        
        dg_values = args[0]
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
                eq = populations[i+1] - populations[i] * np.exp(dg_values[i] / RT) * free_ligand
            elif model == 'nme':
                # Custom model (e.g., n4 depends on n3)
                # ... (extend for other models)
                pass
            eqs.append(eq)
        
        return [eq_metal, eq_ligand] + eqs
    
    # no dg_std provided
    if dg_std == None:
        solution = fsolve(equations, initial_guess, args=dg_values)
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
    # with dg_std: Monte Carlo approach
    elif type(dg_std) == list and len(dg_std)==n_steps:        
        num_samples = 10000
        logk_samples = np.zeros((num_samples, n_steps))
        np.random.seed(42)

        for n in range(num_samples):
            # generate sample values for dg
            dg_samples = np.random.normal(dg_values, dg_std)

            solution = fsolve(equations, initial_guess, args=dg_samples)
            populations = solution[:n_steps + 1]
            free_ligand = solution[-1]
    
            # Compute sample logK values
            for i in range(n_steps):
                if populations[i] == 0:
                    logk_samples[n, i] = np.nan
                else:
                    k = (populations[i+1] / populations[i]) / free_ligand
                    logk_samples[n, i] =  np.log10(k)
            
            # calculate logk mean and std
        logk_mean = np.mean(logk_samples, axis=0)
        logk_std = np.std(logk_samples, axis=0)

        return logk_mean.tolist(), logk_std.tolist()
    else:
        raise TypeError('Std values for ΔG must be a list with the same length as ΔG values.')



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
    result = solve_equilibrium(total_metal, total_ligand, volume_nm3, logk_values)
    fmetal = result['free_metal'] 
    populations = np.append(fmetal,
                             np.array(result['complexes']))
    # transform populations being adimensional
    populations = populations * (0.6022*volume_nm3)

    free_ligand = result['free_ligand'] * (0.6022*volume_nm3)
    
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