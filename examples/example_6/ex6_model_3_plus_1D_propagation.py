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

import h5py
import json
from typing import Dict, Optional

from pycutwm.compile_package import compile_cutwm
from pycutwm.core import run_cutwm_sweep
import pycutwm.plot_tools as pltcr
from pycutwm.power_xy import compute_power_vs_time_from_h5, plot_power_vs_time
import pycutwm.ref_ind_crystals as ri
from pycutwm.gif_from_H5 import create_gif



# Read JSON file
json_path="config.json"
with open(json_path, 'r') as f:
        config = json.load(f)


# Package compilation and execution. 
# Change to "False" if compilation is no required
compile = False
if compile:
    compile_cutwm(dispersion=True,
                  twm_process="opg",
                  cuTWM_path='../../src/',
                  config_path=json_path,
                  compiler="nvcc", opt=3)


# The parameter name in the "param" variable 
# exactly matches the name in the JSON file. 
param = "pump_power_W"
pump_power = 1000000.0
run_code = False
if run_code:
        run_cutwm_sweep(config_path=json_path,
                        param_path="fields/pump/"+param,
                        values=[pump_power],
                        exe_path="../../src/./cuTWM",
                        out_dir_base="results/pump_power",
                        keep_config=True,
                        source_dir=".",             
                        value_decimals=1
        )


# Take values from JSON file
lp = config["crystal"]["wavelengths"]["pump"]
ls = config["crystal"]["wavelengths"]["signal"]
li = config["crystal"]["wavelengths"]["idler"]
temp = 27.0
n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)
n_dict={'pump': n_p, 'signal': n_s, 'idler': n_i}


# --- Example usage ---
field = "pump"  # "pump", "signal", or "idler"
path_to_plot = f"results/pump_power/simulation_pump_power_W_{pump_power}/"
file_to_plot = f"{field}_input_XY.h5"
times, powers = compute_power_vs_time_from_h5(
    h5_path= path_to_plot+file_to_plot,
    time_h5_path= path_to_plot+"time.h5",
    field= field,
    n_dict=n_dict,
    config_path= json_path
)

plot_power_vs_time(times, powers, field=field, figname=f"input_pump_vs_time")


file_to_plot = f"{field}_output_XY.h5"
times, powers = compute_power_vs_time_from_h5(
    h5_path= path_to_plot+file_to_plot,
    time_h5_path= path_to_plot+"time.h5",
    field= field,
    n_dict=n_dict,
    config_path= json_path
)

plot_power_vs_time(times, powers, field=field, figname=f"output_pump_vs_time")


field = "signal"  # "pump", "signal", or "idler"
file_to_plot = f"{field}_output_XY.h5"
times, powers = compute_power_vs_time_from_h5(
    h5_path= path_to_plot+file_to_plot,
    time_h5_path= path_to_plot+"time.h5",
    field= field,
    n_dict=n_dict,
    config_path= json_path
)

plot_power_vs_time(times, powers, field=field, figname=f"output_signal_vs_time")


field = "idler"  # "pump", "signal", or "idler"
file_to_plot = f"{field}_output_XY.h5"
times, powers = compute_power_vs_time_from_h5(
    h5_path= path_to_plot+file_to_plot,
    time_h5_path= path_to_plot+"time.h5",
    field= field,
    n_dict=n_dict,
    config_path= json_path
)

plot_power_vs_time(times, powers, field=field, figname=f"output_idler_vs_time")


## Create GIFs of the field propagation



# === SetupConfiguración ===
field   = "signal"  # field to visualize: "pump", "signal" o "idler"
h5_path = f"results/pump_power/simulation_pump_power_W_1000000.0/{field}_output_XY.h5"   # <- pon acá tu archivo HDF5
gif_name = f"{field}_animation_3D_from_h5.gif"
fps = 10  # frames per second

create_gif(h5_path=h5_path, gif_name=gif_name, field=field, fps=fps, n_dict=n_dict)