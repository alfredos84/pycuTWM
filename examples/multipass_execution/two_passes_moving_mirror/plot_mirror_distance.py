import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("efficiency_vs_distance.dat")
distances_cm = data[:,0]
efficiencies = data[:,1]    

fig = plt.figure(figsize=(8,6))
plt.plot(distances_cm, efficiencies, label='SHG Efficiency', color='blue', markersize=4, marker='o')
plt.xlabel('Mirror Distance (cm)', fontsize=14)
plt.ylabel('SHG Efficiency (%)', fontsize=14)
plt.title('SHG Efficiency vs Mirror Distance', fontsize=16)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.show()
# fig.savefig("efficiency_moving_mirror.png", dpi=150, bbox_inches='tight')
print(np.max(efficiencies))


plt.figure(figsize=(8,6))
plt.plot(distances_cm, efficiencies/np.max(efficiencies), label='SHG Efficiency', color='blue', markersize=4, marker='o')
plt.xlabel('Mirror Distance (cm)', fontsize=14)
plt.ylabel('SHG Efficiency (%)', fontsize=14)
plt.title('SHG Efficiency vs Mirror Distance', fontsize=16)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.show()