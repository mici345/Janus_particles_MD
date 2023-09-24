# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

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

def plot_PMF_same_graph(coord_names, start_frame):
    #######################################
    ### Load data
    #######################################
    folder = "./md/"
    cols = list(range(22))
    start_frame_str = str(start_frame)
    figure_name = folder + 'figures/pmf-plot-' + start_frame_str + '.png'

    # colors for histogram
    colors = cm.get_cmap('viridis', len(coord_names))
    y_min = y_max = 0

    for i in range(0, len(coord_names)):
        coord_name = coord_names[i]
        profile_filename = folder + "profile-" + coord_name + "-" + start_frame_str + ".xvg"

        x_profile, y_profile = read_xvg_file_2cols(profile_filename)

        #######################################
        ### Plot only the PMF
        #######################################
        n_running = 9
        n_shift = int((n_running - 1) / 2)
        y_running_mean = running_mean(y_profile, n_running)
        y_min = min(y_min, min(y_running_mean))
        y_max = max(y_max, max(y_running_mean))
        # plot the graph
        # plt.plot(x_profile, y_profile, color=colors(i))#, alpha=0.3)
        # plot the running average
        plt.plot(x_profile[n_shift:-n_shift], y_running_mean, color=colors(i), label = 'PMMA/' + coord_name)
        plt.xlabel("$\\xi$ ($nm$)") #Distance from PS CoM to PMMA CoM
        plt.ylabel("PMF ($kcal$ $mol^{-1}$)")
    plt.legend()
    plt.title("Potential of Mean Force")

    y_min -= 0.15
    y_max += 0.15
    plt.ylim([y_min, y_max])

    # average surface radius
    plt.plot([1.6, 1.6], [y_min, y_max], color='grey', alpha = 0.3)

    plt.gcf().set_dpi(300)
    plt.savefig(figure_name + '.png', dpi = 300)
    plt.show()


def plot_PMF(coord_name, start_frame):
    #######################################
    ### Load data
    #######################################
    # filenames
    folder = "./md/"
    cols = list(range(22))
    start_frame_str = str(start_frame)
    figure_name = folder + 'figures/pmf-plot-' + coord_name + '-' + start_frame_str + '.png'
    profile_filename = folder + "profile-" + coord_name + "-" + start_frame_str + ".xvg"
    histo_filename = folder + "histo-" + coord_name + "-" + start_frame_str + ".xvg"

    x_profile, y_profile = read_xvg_file_2cols(profile_filename)
    histo_data = np.loadtxt(histo_filename, delimiter='\t', skiprows=17, usecols=cols)

    # find min max for distances
    x_min = min(min(x_profile), min(histo_data[:,0]))
    x_max = max(max(x_profile), max(histo_data[:,0]))
    # add space from each side
    x_min = x_min - 0.02 * (x_max - x_min)
    x_max = x_max + 0.02 * (x_max - x_min)
    x_min = 0.8
    x_max = 2.9
    y_min = -3
    y_max = 4.5

    # colors for histogram
    colors = cm.get_cmap('viridis', histo_data.shape[1] + 5)

    y_min = min(y_profile)
    y_max = max(y_profile)
    deltaG = y_max - y_min
    print("delta G = ", deltaG)

    #######################################
    ### Plot the PMF with count histogram
    #######################################

    fig = plt.figure(tight_layout=True)
    gs = gridspec.GridSpec(3, 1)

    ax = fig.add_subplot(gs[0:2, 0])

    ax.set_title("Potential of Mean Force, PS - PMMA/" + coord_name)

    n_running = 9
    n_shift = int((n_running - 1) / 2)
    y_running_mean = running_mean(y_profile, n_running)
    # plot the graph
    ax.plot(x_profile, y_profile, color=colors(3), alpha=0.3)
    # plot the running average
    ax.plot(x_profile[n_shift:-n_shift], y_running_mean, color=colors(3))
    ax.set_ylabel("PMF ($kcal$ $mol^{-1}$)")
    ax.set_xlim([x_min, x_max])
    y_min = -3
    y_max = 4.5
    ax.set_ylim([y_min, y_max])
    # average surface radius
    plt.plot([1.6, 1.6], [y_min, y_max], color='grey', alpha = 0.3)
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    n_running = 3
    n_shift = int((n_running - 1) / 2)
    y_max = 0;
    ax = fig.add_subplot(gs[2, 0])
    for i in range(1, histo_data.shape[1]):
        x = histo_data[:,0]
        y = histo_data[:,i]
        tmp = max(y)
        if tmp > y_max: y_max = tmp
        y_running_mean = running_mean(histo_data[:,i], n_running)
        idx = y > 0
        ax.plot(x[idx], y[idx], color=colors(i), alpha=1)
        # ax.plot(x[n_shift:-n_shift], y_running_mean, color=colors(i), label=i)
    ax.set_xlabel("$\\xi$ ($nm$)") #Distance from PS CoM to PMMA CoM
    ax.set_ylabel("Count")
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([0, 1.05 * y_max])
    # ax.legend()

    plt.gcf().set_dpi(300)
    plt.savefig(figure_name + '.png', dpi = 300)
    plt.show()


coord_names = ["C2", "C3", "C15", "C28", "C29"]
start_frame = 750

for coord_name in coord_names: plot_PMF(coord_name, start_frame)

plot_PMF_same_graph(coord_names, start_frame)
