"""
ex5_model_3_plus_1D_propagation.py
------------

This code uses the full model, including scattering and diffraction.
The code outputs the time-domain electric fields.
A focused Gaussian pulse with a FWHM duration of 200 fs and 
a beam waist of 30 μm is used as the pump.

"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import imageio

from pycutwm.core import compile_cutwm, run_cutwm_sweep
import pycutwm.plot_tools as pltcr
import pycutwm.ref_ind_crystals as ri


# Read JSON file
json_path="config.json"
with open(json_path, 'r') as f:
        config = json.load(f)


# Package compilation and execution. 
# Change to "False" if compilation is no required
compile = True
if compile:
    compile_cutwm(dispersion=True,
                  twm_process="opg",
                  cuTWM_path='../../src/',
                  config_path=json_path,
                  compiler="nvcc", opt=3)


# The parameter name in the "param" variable 
# exactly matches the name in the JSON file. 
param = "pump_power_W"
pump_power = 10000000.0
run_cutwm_sweep(config_path=json_path,
                param_path="fields/pump/"+param,
                values=[pump_power],
                exe_path="../../src/./cuTWM",
                out_dir_base="results/pump_power",
                keep_config=True,
                source_dir=".",             
                value_decimals=1
)


# # Take values from JSON file
# lp = config["crystal"]["wavelengths"]["pump"]
# ls = config["crystal"]["wavelengths"]["signal"]
# li = config["crystal"]["wavelengths"]["idler"]
# temp = 27.0
# n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)
# n_dict={'pump': n_p, 'signal': n_s, 'idler': n_i}



# # Plot GIF of beam profiles in plane x-y inside the crystal  
# pltcr.field_prop_in_time_make_gif(
#     files_folder="results/pump_power/simulation_"+param+"_"+str(pump_power)+"/", # Here pump power is 10 W
#     config_path=json_path,
#     n_dict=n_dict,
#     gif_path="gifs/propag.gif",
#     cmap='inferno',
#     fps=10,
#     vmin=None,
#     vmax=None
# )