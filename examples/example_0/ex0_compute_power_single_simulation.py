"""
ex0_compute_power_single_simulation.py
------------

This example calculates the output power for a single simulation.

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
    compile_cutwm(dispersion=False, twm_process="sfg",
        cuTWM_path='../../src/', config_path=json_path, 
        compiler="nvcc", opt=3)


# The parameter name in the "param" variable 
# exactly matches the name in the JSON file. 
param = "pump_power_W"
pump_powers = [2.0]
run_code = True
if run_code:
    run_cutwm_sweep(config_path=json_path,
                    param_path="fields/pump/"+param,
                    values=pump_powers,
                    exe_path="../../src/./cuTWM",
                    out_dir_base="results/pump_power",
                    keep_config=True,
                    source_dir=".",             
                    value_decimals=1
    )

# Take values from JSON file
in_power =  config["fields"]["pump"]["pump_power_W"]
lp = config["crystal"]["wavelengths"]["pump"]
ls = config["crystal"]["wavelengths"]["signal"]
li = config["crystal"]["wavelengths"]["idler"]
temp = 27.0
n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)

n_dict={'pump': n_p, 'signal': n_s, 'idler': n_i}

results = []
p_powers,s_powers = [], []

for pump in pump_powers:
    # Build folder name with consistent formatting (same as your sweep output)
    folder = f"results/pump_power/simulation_pump_power_W_{pump:.1f}"
    print(folder)
    # Compute output powers
    powers = compute_power_from_h5(folder=folder, field="p", n=n_p, config_path=json_path)
    p_powers.append(powers["power"])
    print(f"Pump input: {pump:.1f} W, output powers: {powers}")
    powers = compute_power_from_h5(folder=folder, field="s", n=n_s, config_path=json_path)
    s_powers.append(powers["power"])
    results.append({"pump_input": pump, **powers})
    print("Efficiency = "+str(100*s_powers[-1]/in_power))