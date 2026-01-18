import os
import json
import numpy as np
import matplotlib.pyplot as plt
from pycutwm.power_xy import compute_power_from_h5
import pycutwm.ref_ind_crystals as ri

# --- Parámetros globales ---
base_dir = "results_Lcr"
Lcr_values = [10, 20, 30, 40, 50]   # mm
n_passes = 6
pump_powers = [0.104, 0.205, 0.405, 0.6, 0.803, 1.0, 1.207, 1.401, 1.604,
               1.802, 2.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000]

config_path = "config.json"

# --- Índice de refracción (solo una vez) ---
with open(config_path, 'r') as f:
    config = json.load(f)
ls = config["crystal"]["wavelengths"]["signal"]
temp = config["crystal"]["properties_pp"]["temperature"]
n_s = ri.n_mgoppln(ls, temp)
sim_folder = "/home/alfredo/pycuTWM/examples/example_8/results_Lcr/results_Lcr_20mm/pump_power_2_passes/simulation_pump_power_W_2.000"
powers = compute_power_from_h5(folder=sim_folder, field="s", n=n_s, config_path=config_path)
print(powers)