"""
test_performance.py
------------

This code measures the performance for different grid sizes
according to the installed GPU

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


# Read JSON file
json_path="config.json"
with open(json_path, 'r') as f:
        config = json.load(f)


nx_list = [64, 128]
nz_list = [100000]
nt_list = [256, 512, 1024, 2048, 4096]
# nx_list = [64]
# nz_list = [1000]
# nt_list = [256]

for nz in nz_list:
    config["grid"]["grid_points"]["NZ"] = nz
    with open(json_path, 'w') as f:
        json.dump(config, f, indent=4)
    for nt in nt_list:
        config["grid"]["grid_points"]["NT"] = nt
        with open(json_path, 'w') as f:
            json.dump(config, f, indent=4)
        for nx in nx_list: 
            config["grid"]["grid_points"]["NX"] = nx
            config["grid"]["grid_points"]["NY"] = nx
            with open(json_path, 'w') as f:
                json.dump(config, f, indent=4)

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
            pump_power = 1000000.0
            run_code = True
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
