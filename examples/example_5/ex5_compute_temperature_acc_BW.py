"""
ex4_compute_efficiency_mismatch.py
------------

This example calculates the normalized SHG efficiency as a function
of the mismatch factor, Δk. When working with focused Gaussian beams,
maximum efficiency is usually not reached when Δk=0 due to the Gouy phase.
When calibrating an experiment to the maximum output intensity, 
it is important to know the value of Δk in the simulations.

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
compile=False
if compile:
    compile_cutwm(dispersion=False, twm_process="shg",
            cuTWM_path='../../src/', config_path=json_path, 
            compiler="nvcc", opt=3)

# The parameter name in the "param" variable 
# exactly matches the name in the JSON file. 
param = "dk"
temperature = np.arange(52, 64, 0.1)  # in Celsius
dk = []
for temp in list(temperature):  # Define dk for each temperature
    # Take values from JSON file
    lp = config["crystal"]["wavelengths"]["pump"]
    ls = config["crystal"]["wavelengths"]["signal"]
    li = config["crystal"]["wavelengths"]["idler"]
    Lambda = config["crystal"]["properties_pp"]["grating_period"]
    n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)
    dk.append(4*np.pi/lp*(n_s - n_p)- 2*np.pi/Lambda)  # in rad/mm


run_code = False
if run_code:
    run_cutwm_sweep(config_path=json_path,
                    param_path="crystal/wavelengths/"+param,
                    values=dk,
                    exe_path="../../src/./cuTWM",
                    out_dir_base="results/mismatch",
                    keep_config=True,
                    source_dir=".",             
                    value_decimals=7
    )


# Take values from JSON file
waist_um =  config["fields"]["pump"]["waist_um"]
in_power =  config["fields"]["pump"]["pump_power_W"]
Lcr = config["crystal"]["dimensions"]["Lcr"]
lp = config["crystal"]["wavelengths"]["pump"]
ls = config["crystal"]["wavelengths"]["signal"]
li = config["crystal"]["wavelengths"]["idler"]
temp = 27.0
n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)


results = []
p_powers,s_powers = [], []

for mismatch in list(dk):
    # Build folder name with consistent formatting (same as your sweep output)
    folder = f"results/mismatch/simulation_dk_{mismatch:.7f}"
    # Compute output powers
    powers = compute_power_from_h5(folder=folder, field="s", n=n_s, config_path=json_path)
    s_powers.append(powers["power"])

    # results.append({"pump_input": mismatch, **powers})

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

output = np.array(s_powers)

fig, ax = plt.subplots(figsize=(5, 5))

x, y = temperature, 100*output/in_power
ax.plot(x, y/np.max(y), marker='o', markersize=3, linewidth=1, label='SHG efficiency')

ax.set_xlabel(r"Temperature ($^{\circ}$C)")
ax.set_ylabel(r"Normalized SHG efficiency (a.u.)")
ax.grid(True, which='both', linestyle='--', alpha=0.6)
# ax.legend(loc='best')
plt.tight_layout()
plt.show()
fig.savefig("temperature_acc_BW.png", dpi=150, bbox_inches='tight')
plt.close(fig)
np.savetxt('temperature_acc_BW.dat', np.column_stack((x, y)), fmt='%.6f')

# imax = np.argmax(y)
# x_max, y_max = x[imax],y[imax]

# print("Max efficiency at dk = "+ str(x_max))