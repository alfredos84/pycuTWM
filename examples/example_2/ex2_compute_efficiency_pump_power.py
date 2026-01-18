"""
ex2_compute_efficiency_pump_power.py
------------

This example calculates the power scaling for the SHG process.

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


# Datos experimentales (pases 1-5) en W (última columna añadida, convertido de mW a W)
# Borrar esto y la lista de potencias experimentales para volver al ejemplo anterior.

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
2 0,1416 0,453 0,742 0,936 0.844
"""
exp_arr = np.array([[float(num.replace(',', '.')) for num in line.split()] for line in exp_data.strip().split('\n')])

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

for npass in range(1,2):  # From 2 to 5 passes
    # The parameter name in the "param" variable 
    # exactly matches the name in the JSON file. 
    param = "pump_power_W"
    pump_powers = range(1,31)  # From 1 W to 30 W
    config["mul_pass_scheme"]["npasses"] = npass
    if (npass == 1):
        config["mul_pass_scheme"]["multipass"] = False
    else:
        config["mul_pass_scheme"]["multipass"] = True
    with open(json_path, 'w') as f:
        json.dump(config, f, indent=4)

    run_code = False
    if run_code:
        run_cutwm_sweep(config_path=json_path,
                        param_path="fields/pump/"+param,
                        values=pump_powers,
                        exe_path="../../src/./cuTWM",
                        out_dir_base=f"results_{npass}_passes/pump_power",
                        keep_config=True,
                        source_dir=".",             
                        value_decimals=3
        )

    # Take values from JSON file
    lp = config["crystal"]["wavelengths"]["pump"]
    ls = config["crystal"]["wavelengths"]["signal"]
    li = config["crystal"]["wavelengths"]["idler"]
    temp = 27.0
    n_p, n_s, n_i = ri.n_mgoppln(lp,temp), ri.n_mgoppln(ls,temp), ri.n_mgoppln(li,temp)


    n_dict={'pump': n_p, 'signal': n_s, 'idler': n_i}

    results = []
    p_powers,s_powers = [], []

    for pump in pump_powers:
        # Build folder name with consistent formatting (same as your sweep output)
        folder = f"results_{npass}_passes/pump_power/simulation_pump_power_W_{pump:.3f}"
        print(folder)
        # Compute output powers
        powers = compute_power_from_h5(folder=folder, field="p", n=n_p, config_path=json_path)
        p_powers.append(powers["power"])
        print(f"Pump input: {pump:.1f} W, output powers: {powers}")
        powers = compute_power_from_h5(folder=folder, field="s", n=n_s, config_path=json_path)
        s_powers.append(powers["power"])
        results.append({"pump_input": pump, **powers})



    ###################### PLOT ######################

    # Use LaTeX fonts for all
    plt.rcParams.update({
        "font.family": "serif",
        "axes.labelsize": 16,
        "font.size": 14,
        "axes.titlesize": 14,
        "legend.fontsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })

    efficiency = 100 * np.array(s_powers) / np.array(pump_powers)
    pump_powers= np.array(pump_powers)
    s_powers = np.array(s_powers)

    filename = f"Power_scaling_{npass}_passes.txt"
    print("Guardando archivo:", filename)
    np.savetxt(filename, np.column_stack((pump_powers, s_powers)), fmt='%.6f')
    
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(pump_powers, efficiency, marker='o', linewidth=2, label='SHG efficiency')
    # efficiency_exp = exp_arr[:, npass]*100/exp_arr[:, 0]
    # ax.plot(exp_arr[:, 0], efficiency_exp, 'x', color='blue', ms=11, mew=2.3, label='Experiment')
    ax.set_xlabel(r"Pump power (W)")
    ax.set_ylabel(r"Efficiency, $\eta$ (%)")
    ax.set_title(r"SHG efficiency for MgO:PPLN crystal")
    ax.grid(True, which='both', linestyle='--', alpha=0.6)
    # ax.legend(loc='best')
    ax.set_ylim(0, 50)
    plt.tight_layout()
    # plt.show()
    fig.savefig(f"Power_scaling_{npass}_passes.png", dpi=150, bbox_inches='tight')
    plt.close(fig)

    