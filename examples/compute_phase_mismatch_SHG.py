"""
compute_phase_mismatch_SHG.py
------------
This code is used to calculate and graph the mismatch factor, Δk,
as a function of different parameters. This depends on the type
of nonlinear crystal used. For quasi-phase matched crystals, the
temperature and grating period are used. For birefringent crystals,
the polarization of each beam is used.

"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import imageio

from pycutwm.core import compile_cutwm, run_cutwm_sweep
from pycutwm.power_xy import compute_powers_from_folder
import pycutwm.plot_tools as pltcr
import pycutwm.ref_ind_crystals as ri


json_path="./config.json" # read json file
with open(json_path, 'r') as f:
        config = json.load(f)


Lambda = np.arange(5, 10, 0.1, dtype=float) # config["crystal"]["properties_pp"]["grating_period"]
lp = config["crystal"]["wavelengths"]["pump"]
ls = config["crystal"]["wavelengths"]["signal"]
li = config["crystal"]["wavelengths"]["idler"]
temp = 27 #np.array(range(20,80))

# Compute Δk
dk = 4*np.pi/lp*(ri.n_mgoppln(ls,temp)-ri.n_mgoppln(lp,temp))-2*np.pi/Lambda

########### PLOT ###########
# Use LaTeX fonts for all
plt.rcParams.update({
    "font.family": "serif",
    "axes.labelsize": 16,
    "font.size": 14,
    "axes.titlesize": 18,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(Lambda, np.abs(dk), marker='o', linewidth=2, label='$\Delta k$')
print(Lambda[np.argmin(np.abs(dk))])
ax.set_xlabel(r"Tempretature ($^{\circ}$C)")
ax.set_ylabel(r"$\Delta k$ ($\mu$m$^{-1}$)")
ax.grid(True, which='both', linestyle='--', alpha=0.6)
ax.legend(loc='best')

plt.tight_layout()
plt.show()