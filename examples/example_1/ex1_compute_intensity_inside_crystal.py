"""
ex1_compute_intensity_inside_crystal.py
------------

This code calculates the electric fields inside the nonlinear crystal (not just at the output).
It can be useful for calculating the intensity in units of W/cm² and studying the damage threshold
when working with focused beams. In the JSON file, you can vary the focal position or other
parameters and study the generation of the fields. The intensity profiles are plotted in
a cut in the y-z plane.

"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import imageio

from pycutwm.compile_package import compile_cutwm
from pycutwm.core import run_cutwm_sweep
import pycutwm.plot_tools as pltcr
import pycutwm.ref_ind_crystals as ri


# Read JSON file
json_path="config.json"
with open(json_path, 'r') as f:
        config = json.load(f)


# Package compilation and execution. 
# Change to "False" if compilation is no required
compile = False
if compile:
    compile_cutwm(dispersion=False, twm_process="shg",
        cuTWM_path='../../src/', config_path=json_path, 
        compiler="nvcc", opt=3)


# The parameter name in the "param" variable 
# exactly matches the name in the JSON file. 
param = "pump_power_W"
pump_power = [10.0]
run_code = False
if run_code:
        run_cutwm_sweep(config_path=json_path,
                        param_path="fields/pump/"+param,
                        values=pump_power,
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


pltcr.plot_intensity_cut_YZ_from_h5(
    files_folder="results/pump_power/simulation_"+param+"_10.0/", # Here pump power is 10 W
    config_path=json_path,
    n_dict=n_dict,
    cmap='inferno',
    image_dir="intensity_cut_YZ/",
    image_name="intensity_profile",
    image_extension=".png",
    y_limits=(-60, 60),     # e.g. (-500, 500)
    z_limits=None      # e.g. (0, 20000)
)