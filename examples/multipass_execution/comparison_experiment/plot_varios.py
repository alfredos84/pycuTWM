import os
import json
import numpy as np
import matplotlib.pyplot as plt
from pycutwm.power_xy import compute_power_from_h5
import pycutwm.ref_ind_crystals as ri

# --- Parámetros globales ---
base_dir = "results_Lcr"
Lcr_values = [20] #[10, 20, 30, 40, 50]   # mm
n_passes = 1
pump_powers = [0.104, 0.205, 0.405, 0.6, 0.803, 1.0, 1.207, 1.401, 1.604,
               1.802, 2.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000]

config_path = "config.json"

# --- Índice de refracción (solo una vez) ---
with open(config_path, 'r') as f:
    config = json.load(f)
ls = config["crystal"]["wavelengths"]["signal"]
temp = config["crystal"]["properties_pp"]["temperature"]
n_s = ri.n_mgoppln(ls, temp)

# --- Función auxiliar ---
def find_sim_folder(base_folder, p):
    """Devuelve la ruta de la carpeta de simulación para una potencia p."""
    folder = os.path.join(base_folder, f"simulation_pump_power_W_{p:.3f}")
    if not os.path.exists(folder):
        # intentar sin formato fijo (por si no tiene ceros a la derecha)
        folder = os.path.join(base_folder, f"simulation_pump_power_W_{p}")
    return folder if os.path.exists(folder) else None

# --- Loop principal ---
for Lcr in Lcr_values:
    Lcr_folder = os.path.join(base_dir, f"results_Lcr_{Lcr}mm")
    if not os.path.exists(Lcr_folder):
        print(f"[WARN] Carpeta no encontrada: {Lcr_folder}")
        continue

    print(f"\nProcesando Lcr = {Lcr} mm ...")

    plt.figure(figsize=(9,6))
    colors = plt.cm.viridis(np.linspace(0, 1, n_passes))

    for pass_idx in range(1, n_passes+1):
        pass_folder = os.path.join(Lcr_folder, f"pump_power_{pass_idx}_passes")
        if not os.path.exists(pass_folder):
            print(f"  [WARN] Falta {pass_folder}")
            continue

        y_signal = []
        eff = []

        for p in pump_powers:
            sim_folder = find_sim_folder(pass_folder, p)
            if sim_folder is None:
                print(f"    [MISS] pump={p} W (no existe carpeta)")
                y_signal.append(np.nan)
                continue

            try:
                powers = compute_power_from_h5(folder=sim_folder, field="s", n=n_s, config_path=config_path)
                P_s = powers["power"]
                eta = P_s / p
                eff.append(eta)
                y_signal.append(P_s)
            except Exception as e:
                print(f"    [ERR] {sim_folder}: {e}")
                y_signal.append(np.nan)
                eff.append(np.nan)

        # Graficar solo si hay datos válidos
        eff = np.array(eff)
        if np.isfinite(eff).any():
            plt.plot(pump_powers, eff*100, '-o', color=colors[pass_idx-1],
                     lw=2, label=f"{pass_idx} pass{'es' if pass_idx>1 else ''}")

    # --- Personalización del gráfico ---
    plt.xlabel("Input Pump Power (W)", fontsize=13)
    plt.ylabel("Conversion Efficiency", fontsize=13)
    plt.title(f"SHG Conversion Efficiency vs Pump Power\nLcr = {Lcr} mm", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()

    outname = f"efficiency_Lcr_{Lcr}mm.png"
    plt.savefig(outname, dpi=150)
    print(f"[OK] Figura guardada: {outname}")

plt.show()
