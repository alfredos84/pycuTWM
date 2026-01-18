import os
import numpy as np
import matplotlib.pyplot as plt
import json

from pycutwm.power_xy import compute_power_from_h5
import pycutwm.ref_ind_crystals as ri

# --- Parámetros ---
n_passes = 5
pump_powers = [0.104, 0.205, 0.405, 0.6, 0.803, 1.0, 1.207, 1.401, 1.604, 1.802, 2.0]
main_folders = [f"results_passes_{i+1}/pump_power" for i in range(n_passes)]

# Datos experimentales (pases 1-5) en W (última columna añadida, convertido de mW a W)
exp_data = """
0,104 0,000394 0,001464 0,00302 0,00489 0.00507
0,205 0,001582 0,00572 0,01159 0,01871 0.01905
0,405 0,00618 0,02173 0,0437 0,0672 0.0668
0,6 0,01345 0,047 0,0914 0,1369 0.1347
0,803 0,02392 0,0831 0,1561 0,2277 0.2196
1 0,0368 0,1264 0,2312 0,33 0.310
1,207 0,0533 0,1788 0,324 0,447 0.410
1,401 0,0712 0,2381 0,419 0.56 0.512
1,604 0,0926 0,303 0,522 0,681 0.603
1,802 0,1151 0,374 0,629 0,814 0.697
2 0,1416 0,453 0,742 0,936 0.780
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
            print(f"NO SE ENCONTRÓ carpeta para pump_power={p} en {folder}")
            y.append(np.nan)
            continue
        try:
            powers = compute_power_from_h5(folder=sim_folder, field="s", n=n_s, config_path=config_path)
            y.append(1000 * powers["power"])  # de W a mW
        except Exception as e:
            print(f"Error leyendo {sim_folder}: {e}")
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

# --- Guardar archivos .dat ---

# SIMULACIONES: [Input, Sim1, Sim2, ..., Sim5] (en W)
sim_powers_w = sim_signal_powers.T / 1000  # shape (11, 5), de mW a W
sim_table = np.column_stack([pump_powers, sim_powers_w])
np.savetxt("simulation_powers.dat", sim_table, header="InputPumpPower_W Sim1_W Sim2_W Sim3_W Sim4_W Sim5_W", fmt="%.6f")

# EXPERIMENTOS: [Input, Exp1, Exp2, ..., Exp5] (en W)
exp_table = exp_arr  # ya está en W
np.savetxt("experiment_powers.dat", exp_table, header="InputPumpPower_W Exp1_W Exp2_W Exp3_W Exp4_W Exp5_W", fmt="%.6f")

# EFFICIENCY (simulacion): 100*output_signal/input_pump
sim_eff = 100 * sim_powers_w / np.array(pump_powers)[:, None]  # shape (11,5)
sim_eff_table = np.column_stack([pump_powers, sim_eff])
np.savetxt("simulation_efficiency.dat", sim_eff_table, header="InputPumpPower_W Sim1_% Sim2_% Sim3_% Sim4_% Sim5_%", fmt="%.6f")

# EFFICIENCY (experimentos): igual
exp_eff = 100 * exp_table[:, 1:] / exp_table[:, 0:1]  # shape (11,5)
exp_eff_table = np.column_stack([pump_powers, exp_eff])
np.savetxt("experiment_efficiency.dat", exp_eff_table, header="InputPumpPower_W Exp1_% Exp2_% Exp3_% Exp4_% Exp5_%", fmt="%.6f")
