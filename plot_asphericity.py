# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


data_mda = np.loadtxt('./Asphericity_Rgyr_300ns.csv', delimiter=',', skiprows=1)
time_mda = data_mda[:,1] / 1000.0 # from ps to ns

legend = ['PS', 'PS+SDS+MMA+PMMA ($r_c = 3\AA$)']
colors = ['r', 'b', 'g', 'grey']

### Asphericity
n_running = 15
n_shift = int((n_running - 1) / 2)
for i in range(0,2):
    plt.plot(time_mda, data_mda[:,2+i], color=colors[i], alpha=0.3)
    y_running_mean = running_mean(data_mda[:,2+i], n_running)
    plt.plot(time_mda[n_shift:-n_shift], y_running_mean, color=colors[i], label = legend[i])
plt.ylim([0.0, 0.2])
plt.xlabel("Time, ns")
plt.ylabel("Asphericity")
plt.legend()
plt.show()

### Radius of gyration
n_running = 15
n_shift = int((n_running - 1) / 2)
for i in range(0,2):
    plt.plot(time_mda, data_mda[:,4+i], color=colors[i], alpha=0.3)
    y_running_mean = running_mean(data_mda[:,4+i], n_running)
    plt.plot(time_mda[n_shift:-n_shift], y_running_mean, color=colors[i], label = legend[i])
plt.xlabel("Time, ns")
plt.ylabel("Radius of gyration, nm")
plt.legend()
plt.show()
