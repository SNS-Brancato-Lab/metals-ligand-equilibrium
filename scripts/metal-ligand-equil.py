import os.path as osp
import sys

SCRIPT_DIR = osp.abspath(osp.dirname(__file__))
sys.path.insert(0, osp.join(SCRIPT_DIR, "../"))

from thermodynamic.conversions import dg_to_logk, logk_to_dg
from equilibrium.conc_solver import solve_equilibrium
from formatting.utils import print_results, print_equilibrium


# Find first equilibrium concentrations of the system
logk_values = [5.40, 4.47, 2.23]
results = solve_equilibrium(
    moles_metal=10, 
    moles_ligand=30, 
    volume_nm3=302, 
    logk_values=logk_values
)


# ΔG to logK (sequential model)
# TODO add different possibilities to sequential
dg_values = [28.5, 22.1, 12.3]  # Hypothetical ΔG in kJ/mol
logk = dg_to_logk(
    dg_values=dg_values,
    total_metal=10,
    total_ligand=30,
    volume_nm3=302,
    temp=300,
    model="sequential",
)

# logK to ΔG

dg = logk_to_dg(
    logk_values=logk_values, total_metal=10, 
    total_ligand=30, 
    volume_nm3=302, 
    temp=300
)

print_equilibrium(results)
print_results("Computed logK Values", logk, "Step", "")
print_results("Computed ΔG Values", dg, "Step", "kJ/mol")