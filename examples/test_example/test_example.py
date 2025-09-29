"""
test_example.py
------------

This example serves as a first test after installing the package

"""

import os
import json
from pycutwm.compile_package import compile_cutwm

# Read JSON file
json_path="config.json"
with open(json_path, 'r') as f:
        config = json.load(f)


# Package compilation and execution. 
compile_cutwm(dispersion=False, twm_process="shg",
    cuTWM_path='../../src/', config_path=json_path, 
    compiler="nvcc", opt=3)