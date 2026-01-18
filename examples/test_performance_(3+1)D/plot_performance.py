import numpy as np
import matplotlib.pyplot as plt

# Cargar datos: NX NY NT NZ time
data = np.loadtxt("timings_grid_GPU_A100.txt")

NX = data[:, 0].astype(int)
NY = data[:, 1].astype(int)   # igual a NX, pero lo dejamos por claridad
NZ = data[:, 2].astype(int)
NT = data[:, 3].astype(int)
t  = data[:, 4]

Ns  = np.unique(NX)
NZs = np.unique(NZ)

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 12
})

for nz in NZs:
    plt.figure()

    for N in Ns:
        mask = (NX == N) & (NZ == nz)
        if np.any(mask):
            NT_vals = NT[mask]
            t_vals  = t[mask]

            # ordenar por NT para que la curva salga bien
            order = np.argsort(NT_vals)
            NT_vals = NT_vals[order]
            t_vals  = t_vals[order]

            plt.plot(NT_vals, np.log(t_vals), marker="o", label=f"Nx=Ny={N}")

    plt.xlabel(r"$N_T$")
    plt.ylabel("Runtime (s)")
    plt.title(f"Runtime vs $N_T$ for $N_Z = {nz}$")
    plt.xscale("log", base=2)
    # plt.yscale()  # opcional pero muy útil para ver el scaling
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"Performance_NZ_{nz}.png", dpi=150, bbox_inches="tight")
    plt.close()    
