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
from pycutwm.power_xy import compute_power_from_h5
import pycutwm.ref_ind_crystals as ri



# --- Parámetros ---
n_passes = 1
pump_powers = [0.104, 0.205, 0.405, 0.6, 0.803, 1.0, 1.207, 1.401, 1.604, 1.802, 2.000, 3.000, 4.000, 5.000]
# main_folders = [f"results_passes_{i+1}/pump_power" for i in range(n_passes)]
main_folders = [f"pump_power_{i+1}_passes" for i in range(n_passes)]

# Datos experimentales (pases 1-5) en W (última columna añadida, convertido de mW a W)
exp_data = """
0,104 0,000394 0,001464 0,00302 0,00489 0.00522
0,205 0,001582 0,00572 0,01159 0,01871 0.019595
0,405 0,00618 0,02173 0,0437 0,0672 0.0688
0,6 0,01345 0,047 0,0914 0,1369 0.1411
0,803 0,02392 0,0831 0,1561 0,2277 0.2274
1 0,0368 0,1264 0,2312 0,33 0.323
1,207 0,0533 0,1788 0,324 0,447 0.430
1,401 0,0712 0,2381 0,419 0.56 0.530
1,604 0,0926 0,303 0,522 0,681 0.636
1,802 0,1151 0,374 0,629 0,814 0.737
2 0,1416 0,453 0,742 0,936 0.884
"""
exp_arr = np.array([[float(num.replace(',', '.')) for num in line.split()] for line in exp_data.strip().split('\n')])

# --- Índice de refracción solo una vez ---
config_path = "config.json"
with open(config_path, 'r') as f:
    config = json.load(f)
ls = config["crystal"]["wavelengths"]["signal"]
temp = config["crystal"]["properties_pp"]["temperature"]
n_s = ri.n_mgoppln(ls, temp)

def find_sim_folder(base_folder, p):
    sub = f"simulation_pump_power_W_{p:.3f}"
    sim_folder = os.path.join(base_folder, sub)
    if os.path.exists(sim_folder):
        return sim_folder
    return None

# --- Leer las potencias simuladas ---
sim_signal_powers = []
for pass_idx, folder in enumerate(main_folders):
    y = []
    for p in pump_powers:
        sim_folder = find_sim_folder(folder, p)
        if sim_folder is None or not os.path.exists(sim_folder):
            print(f"Folder not found for pump_power={p} en {folder}")
            y.append(np.nan)
            continue
        try:
            powers = compute_power_from_h5(folder=sim_folder, field="s", n=n_s, config_path=config_path)
            y.append(1000 * powers["power"])  # de W a mW
        except Exception as e:
            print(f"Error reading {sim_folder}: {e}")
            y.append(np.nan)
    sim_signal_powers.append(y)
sim_signal_powers = np.array(sim_signal_powers)

# --- Graficar ---
colors = plt.cm.viridis(np.linspace(0, 1, n_passes))
plt.figure(figsize=(9, 6))

for i in range(n_passes):
    label_sim = f"Sim {i+1} pass{'es' if i>0 else ''}"
    label_exp = f"Exp. {i+1} pass{'es' if i>0 else ''}"
    # Graficar simulación
    plt.plot(pump_powers, sim_signal_powers[i], '-', color=colors[i], lw=2, label=label_sim)
    # Graficar experimental (todos los pases)
    plt.plot(exp_arr[:, 0], exp_arr[:, i+1]*1000, 'x', color=colors[i], ms=11, mew=2.3, label=label_exp)

plt.xlabel("Input Pump Power (W)", fontsize=14)
plt.ylabel("Signal Output Power (mW)", fontsize=14)
plt.title("Signal Output vs Pump Power\nSimulation (lines) & Experiment (crosses)", fontsize=14)
plt.grid(True, which='both', alpha=0.2)
plt.legend(fontsize=11, ncol=2)
plt.tight_layout()
plt.show()
