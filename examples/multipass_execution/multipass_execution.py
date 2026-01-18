"""
ex6_multipass_execution.py
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
compile = True
if compile:
    compile_cutwm(dispersion=False, twm_process="shg",
        cuTWM_path='../../src/', config_path=json_path, 
        compiler="nvcc", opt=3)


for npasses in range(1,2):  # From 1 to 6 passes
    # The parameter name in the "param" variable 
    # exactly matches the name in the JSON file. 
    param = "pump_power_W"
    pump_power = [0.104, 0.205, 0.405, 0.6, 0.803, 1.0, 1.207, 1.401, 1.604, 1.802, 2.0, 3.0, 4.0, 5.0]
    config["mul_pass_scheme"]["npasses"] = npasses
    with open(json_path, 'w') as f:
        json.dump(config, f, indent=4)
    run_code = True
    if run_code:
        run_cutwm_sweep(config_path=json_path,
                        param_path="fields/pump/"+param,
                        values=pump_power,
                        exe_path="../../src/./cuTWM",
                        out_dir_base=f"pump_power_{npasses}_passes",
                        keep_config=True,
                        source_dir=".",             
                        value_decimals=3
        )


# # Take values from JSON file
# lp = config["crystal"]["wavelengths"]["pump"]
# ls = config["crystal"]["wavelengths"]["signal"]
# temp = config["crystal"]["properties_pp"]["temperature"]
# n_p, n_s = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp)
# n_dict={'pump': n_p, 'signal': n_s, 'idler': n_s}


# # Plot y-z profile
# pltcr.plot_intensity_cut_YZ_from_h5(
#     files_folder="results/pump_power/simulation_"+param+"_10.0/", # Here pump power is 10 W
#     config_path=json_path,
#     n_dict=n_dict,
#     cmap='inferno',
#     image_dir="intensity_cut_YZ/",
#     image_name="intensity_profile",
#     image_extension=".png",
#     y_limits=(-60, 60),     # e.g. (-500, 500)
#     z_limits=None      # e.g. (0, 20000)
# )