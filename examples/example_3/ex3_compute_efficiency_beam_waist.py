"""
ex3_compute_efficiency_beam_waist.py
------------

This example calculates SHG efficiency by varying the pump's beam waist.
Efficiency is plotted as a function of the focusing parameter, 
according to Boyd-Kleinman theory.

"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import imageio

from pycutwm.compile_package import compile_cutwm
from pycutwm.core import run_cutwm_sweep
from pycutwm.power_xy import compute_power_from_h5
import pycutwm.ref_ind_crystals as ri


# Read JSON file
json_path="config.json"
with open(json_path, 'r') as f:
        config = json.load(f)


# Package compilation and execution. 
# Change to "False" if compilation is no required
compile = True
if compile:
    compile_cutwm(dispersion=False, twm_process="shg",
        cuTWM_path='../../src/', config_path=json_path, 
        compiler="nvcc", opt=3)


# The parameter name in the "param" variable 
# exactly matches the name in the JSON file. 
param = "waist_um"
beam_waists = range(15,80)
run_code = True
if run_code:
    run_cutwm_sweep(config_path=json_path,
                    param_path="fields/pump/"+param,
                    values=beam_waists,
                    exe_path="../../src/./cuTWM",
                    out_dir_base="results/beam_waist",
                    keep_config=True,
                    source_dir=".",             
                    value_decimals=1
    )


# Take values from JSON file
in_power =  config["fields"]["pump"]["pump_power_W"]
Lcr = config["crystal"]["dimensions"]["Lcr"]
lp = config["crystal"]["wavelengths"]["pump"]
ls = config["crystal"]["wavelengths"]["signal"]
li = config["crystal"]["wavelengths"]["idler"]
temp = 27.0
n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)


results = []
p_powers,s_powers = [], []

for waist in beam_waists:
    # Build folder name with consistent formatting (same as your sweep output)
    folder = f"results/beam_waist/simulation_waist_um_{waist:.1f}"
    # Compute output powers
    powers = compute_power_from_h5(folder=folder, field="s", n=n_s, config_path=json_path)
    s_powers.append(powers["power"])
    results.append({"pump_input": waist, **powers})


########### PLOT ###########
# Use LaTeX fonts for all
plt.rcParams.update({
    "font.family": "serif",
    "axes.labelsize": 16,
    "font.size": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

efficiency = 100 * np.array(s_powers) / in_power
b = Lcr*lp/( np.array(beam_waists)**2*2*np.pi*n_p)
fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(b, efficiency, marker='o', linewidth=2, label='SHG efficiency')
ax.set_xlabel(r"Focusing parameter, $b$")
ax.set_ylabel(r"Efficiency, $\eta$ (%)")
ax.set_title(r"SHG efficiency for MgO:PPLN crystal")
ax.grid(True, which='both', linestyle='--', alpha=0.6)
ax.legend(loc='best')
plt.tight_layout()
# plt.show()
plt.savefig("Efficiency.png", dpi=150, bbox_inches="tight")
plt.close(fig)