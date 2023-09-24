# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def read_xvg_file_2cols(xvg_filename):

    #Read data from the xvg file
    with open(xvg_filename,'r') as xvg_file:
        x,y = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]))
                for line in xvg_file 
                if not line.startswith(("#","@")) 
            ])
        )
        #x = x/10
        return [x, y]
    return [[], []]

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


folder = "./"

### AZM heads relative to Water's oxygens
rdf_AZMheads_SOL_xvg = "rdf-AZM_heads-relative-to-Water_O.xvg"
x_AZMheads2SOL, y_AZMheads2SOL = read_xvg_file_2cols(folder + rdf_AZMheads_SOL_xvg)
# plt.plot(x_AZMheads2SOL, y_AZMheads2SOL, color='b', alpha=0.3)
n_running = 7
n_shift = int((n_running - 1) / 2)
y_running_mean = running_mean(y_AZMheads2SOL, n_running)
plt.plot(x_AZMheads2SOL[n_shift:-n_shift], y_running_mean, color='b', label = "SDS heads")
plt.xlim([0, 1])
plt.title("RDF of SDS heads in reference to water")
plt.xlabel("r, nm")
plt.ylabel("g(r)")
plt.show()

rdf_systems = {"SDS heads" : "AZM_heads", "SDS tails" : "AZM_tails", "MMA" : "WVD_heavy", "PMMA" : "6CB_heavy"}
colors = ['r', 'grey', 'b', 'g']

n_running = 49

### RDF in reference to PS center of mass (CoM)

# norm: rdf (default)
i = 0
for system in rdf_systems:
    xvg_filename = "rdf-" + rdf_systems[system] + "-relative-to-PS_CoM.xvg"
    x, y = read_xvg_file_2cols(folder + xvg_filename)
    plt.plot(x, y, color=colors[i], alpha=0.2)
    n_shift = int((n_running - 1) / 2)
    y_running_mean = running_mean(y, n_running)
    plt.plot(x[n_shift:-n_shift], y_running_mean, color=colors[i], label = system)
    i += 1
    
plt.title("Radial distribution in reference to PS center of mass")
plt.xlabel("r, nm")
plt.ylabel("Cumulative number RDF")
plt.legend()
plt.show()

# norm: number_density
i = 0
for system in rdf_systems:
    xvg_filename = "rdf-" + rdf_systems[system] + "-relative-to-PS_CoM-number_density.xvg"
    x, y = read_xvg_file_2cols(folder + xvg_filename)
    plt.plot(x, y, color=colors[i], alpha=0.2)
    n_shift = int((n_running - 1) / 2)
    y_running_mean = running_mean(y, n_running)
    plt.plot(x[n_shift:-n_shift], y_running_mean, color=colors[i], label = system)
    i += 1
    
plt.title("Radial distribution in reference to PS center of mass")
plt.xlabel("r, nm")
plt.ylabel("Number density")
plt.legend()
plt.show()

# no norm
i = 0
for system in rdf_systems:
    xvg_filename = "rdf-" + rdf_systems[system] + "-relative-to-PS_CoM-nonorm.xvg"
    x, y = read_xvg_file_2cols(folder + xvg_filename)
    plt.plot(x, y, color=colors[i], alpha=0.2)
    n_shift = int((n_running - 1) / 2)
    y_running_mean = running_mean(y, n_running)
    plt.plot(x[n_shift:-n_shift], y_running_mean, color=colors[i], label = system)
    i += 1
    
plt.title("Radial distribution in reference to PS center of mass")
plt.xlabel("r, nm")
plt.ylabel("g(r), non-normalized")
plt.legend()
plt.show()

### RDF in reference to PS surface

i = 0
for system in rdf_systems:
    xvg_filename = "rdf-" + rdf_systems[system] + "-relative-to-PS_surf.xvg"
    x, y = read_xvg_file_2cols(folder + xvg_filename)
    plt.plot(x, y, color=colors[i], alpha=0.3)
    n_running = 15
    n_shift = int((n_running - 1) / 2)
    y_running_mean = running_mean(y, n_running)
    plt.plot(x[n_shift:-n_shift], y_running_mean, color=colors[i], label = system)
    i += 1

plt.title("Radial distribution in reference to PS surface")
plt.xlabel("r, nm")
plt.ylabel("g(r), non-normalized")
plt.legend()
plt.show()
