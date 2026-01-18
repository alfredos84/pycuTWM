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

# from pycutwm.compile_package import compile_cutwm
# from pycutwm.core import run_cutwm_sweep
# from pycutwm.power_xy import compute_power_from_h5
import pycutwm.ref_ind_crystals as ri


Lambda = 6.92 #np.arange(5, 10, 0.1, dtype=float) # config["crystal"]["properties_pp"]["grating_period"]
lp = 1.064 # config["crystal"]["wavelengths"]["pump"]
ls = 0.532 # config["crystal"]["wavelengths"]["signal"]
li = ls #config["crystal"]["wavelengths"]["idler"]
temp = np.arange(55, 60, 0.05, dtype=float)

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

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(temp, np.abs(dk), marker='o', linewidth=2, label='$\Delta k$')
print(temp[np.argmin(np.abs(dk))])
ax.set_xlabel(r"Tempretature ($^{\circ}$C)")
ax.set_ylabel(r"$\Delta k$ ($\mu$m$^{-1}$)")
ax.grid(True, which='both', linestyle='--', alpha=0.6)
ax.legend(loc='best')

plt.tight_layout()
plt.show()
fig.savefig("temperature_acc_BW.png", dpi=150, bbox_inches='tight')
plt.close(fig)