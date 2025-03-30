import os.path as osp
import sys


SCRIPT_DIR = osp.abspath(osp.dirname(__file__))
sys.path.insert(0, osp.join(SCRIPT_DIR, "../"))

from thermodynamic.conversions import dg_to_logk, logk_to_dg
from equilibrium.conc_solver import solve_equilibrium




# ΔG to logK (sequential model)
dg_values = [28.5, 22.1, 12.3]  # Hypothetical ΔG in kJ/mol
logk = dg_to_logk(
    dg_values=dg_values,
    total_metal=10,
    total_ligand=30,
    volume_nm3=302,
    temp=300,
    model='sequential'
)
print("Computed logK:", [f"{k:.2f}" for k in logk])

# logK to ΔG
logk_values = [5.4, 4.47, 2.23]

equils = solve_equilibrium(

    moles_metal=10,
    moles_ligand=30,
    volume_nm3=302,
    k_values=logk_values
    
)
print(equils)

dg = logk_to_dg(
    logk_values=logk_values,
    total_metal=10,
    total_ligand=30,
    volume_nm3=302,
    temp=300
)
print("Computed ΔG:", [f"{g:.2f} kJ/mol" for g in dg])