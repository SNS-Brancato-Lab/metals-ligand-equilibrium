import numpy as np
from scipy.optimize import fsolve

def moles_to_concentration(moles, volume_nm3):
    """Convert number of moles (particles) to concentration (Molar)."""
    # 1 M = 1 mole / liter, 1 liter = 1e27 nm³
    # Avogadro's number = 6.022e23 particles/mol
    return moles / (0.6022 * volume_nm3)  # 0.6022 ≈ 6.022e23 / 1e27 * 1e3

def solve_equilibrium(moles_metal, moles_ligand, volume_nm3, logk_values, temp=300, initial_guess=None):
    """
    Solves equilibrium concentrations for metal-ligand complexes.
    
    Args:
        moles_metal (float): Initial moles of metal ions.
        moles_ligand (float): Initial moles of ligand molecules.
        volume_nm3 (float): Volume of the system in nm³.
        k_values (list): Stepwise equilibrium constants [k1, k2, ...].
        temp (float): Temperature in Kelvin.
        initial_guess (list): Optional initial guess for solver.
    
    Returns:
        dict: Equilibrium concentrations of free metal, ligand, and complexes.
    """
    en0 = moles_to_concentration(moles_ligand, volume_nm3)
    cd0 = moles_to_concentration(moles_metal, volume_nm3)
    n_complexes = len(logk_values)
    k_values = np.array([10**i for i in logk_values]) 
    
    # Define the system of equations dynamically
    def equations(vars):
        en = vars[0]
        cd = vars[1]
        complexes = vars[2:2 + n_complexes]
        
        # Ligand conservation: en0 = en + Σ(i * complex_i)
        eq_en = en0 - en - sum((i+1)*c for i, c in enumerate(complexes))
        # Metal conservation: cd0 = cd + Σ(complex_i)
        eq_cd = cd0 - cd - sum(complexes)
        eqs = [eq_en, eq_cd]
        
        # Complex formation equations
        for i in range(n_complexes):
            if i == 0:
                # c1 = k1 * cd * en
                eq = complexes[i] - k_values[i] * cd * en
            else:
                # ci = ki * c_{i-1} * en
                eq = complexes[i] - k_values[i] * complexes[i-1] * en
            eqs.append(eq)
    
        return(eqs)
    
    # Initial guess
    if initial_guess is None:
        initial_guess = [en0, cd0] + [0.01] * n_complexes
    
    # Solve
    solution = fsolve(equations, initial_guess)
    
    # Extract results
    en_eq, cd_eq = solution[0], solution[1]
    complexes_eq = solution[2:]
    return {
        'free_ligand': en_eq, #*(0.6022 * volume_nm3),
        'free_metal': cd_eq, #*(0.6022 * volume_nm3),
        'complexes': complexes_eq #*(0.6022 * volume_nm3)
    }