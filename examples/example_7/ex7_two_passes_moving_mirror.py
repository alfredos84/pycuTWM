"""
ex7_two_passes_moving_mirror.py
------------

This code computes the SHG efficiency in a double-pass scheme with an external concave mirror
with variable crystal-mirror distance. 
In the JSON file, you can vary mirror distance in cm an plot results at the end of this script.

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
param = "distances_in_air"
distances_cm = np.linspace(3.5, 10.5, 500)
values = [[float(x)] for x in distances_cm]
run_code = True
if run_code:
    run_cutwm_sweep(config_path=json_path,
                    param_path="mul_pass_scheme/"+param,
                    values=values,
                    exe_path="../../src/./cuTWM",
                    out_dir_base="results/distances_in_air",
                    keep_config=True,
                    source_dir=".",             
                    value_decimals=3
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

for d in list(distances_cm):
    # Build folder name with consistent formatting (same as your sweep output)
    if isinstance(d, (list, np.ndarray)):
        d_str = "__".join([f"{v:.3f}" for v in d])
    else:
        d_str = f"{d:.3f}"
    folder = f"results/distances_in_air/simulation_distances_in_air_{d_str}"
    
    # Compute output powers
    powers = compute_power_from_h5(folder=folder, field="s", n=n_s, config_path=json_path)
    s_powers.append(powers["power"])
    print("Signal power = "+str(s_powers[-1])+ " W")
    # results.append({"pump_input": mismatch, **powers})

########### PLOT ###########
# Use LaTeX fonts for all
plt.rcParams.update({
    "font.family": "serif",
    "axes.labelsize": 16,
    "font.size": 14,
    "axes.titlesize": 18,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

output = np.array(s_powers)

fig, ax = plt.subplots(figsize=(5, 5))

x, y = distances_cm, 100*output/in_power
ax.plot(x, y, marker='o', linewidth=1, label='SHG efficiency', color='blue')
plt.ylim([0,30])
plt.vlines(x=9.2, ymin=0, ymax=25, colors='black', linestyles='dashed', linewidth=2)
plt.vlines(x=4.3, ymin=0, ymax=25, colors='black', linestyles='dashed', linewidth=2)
ax.set_xlabel(r"Mirror distance (cm)")
ax.set_ylabel(r"SHG efficiency (a.u.)")
ax.set_title(r"SHG efficiency for MgO:PPLN crystal")
ax.grid(True, which='both', linestyle='--', alpha=0.6)
ax.legend(loc='best')
plt.tight_layout()
plt.show()
fig.savefig("efficiency_distance_mirror_waist_"+str(waist_um)+".png", dpi=150, bbox_inches='tight')
plt.close(fig)

np.savetxt('efficiency_vs_distance.dat', np.column_stack((x, y)), fmt='%.6f')