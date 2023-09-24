# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

#########
#PMMA contacts with MMA and SDS and PS:
#- distribution of the number of PMMA atoms within 3 A of PS atoms, and a version normalized by the number of PS atoms,
#- distribution of the number of PMMA atoms within 3 A of MMA atoms, and a version normalized by the number of MMA atoms,
#- distribution of the number of PMMA atoms within 3 A of SDS atoms, and a version normalized by the number of SDS atoms,
#- a plot of the total number of MMA molecules and the total number of MMA atoms that are within 3 A of PMMA over the trajectory to show the evolution of the connections (and a normalized version for the number of atoms),
#- a plot of the total number of SDS molecules and the total number of SDS atoms that are within 3 A of PMMA over the trajectory to show the evolution of the connections (and a normalized version for the number of atoms).
#########

rc = 2.5 # angstroms
data_mda = np.loadtxt('./data/stats_mda_frames_0-1500_'+ str(rc) + 'A.csv', delimiter = ',', skiprows = 1)
rcs = str(rc) + '$\AA$'

headers = [
    'id', 'frame', 'time (ps)',
    'SDS atoms within r_c of PS', 'SDS compounds within r_c of PS',
    'MMA atoms within r_c of PS', 'MMA compounds within r_c of PS',
    'PMMA atoms within r_c of PS', 'PMMA compounds within r_c of PS',
    'PMMA atoms within r_c of surface SDS', 'PMMA compounds within r_c of surface SDS',
    'PMMA atoms within r_c of surface MMA', 'PMMA compounds within r_c of surface MMA',
    'SDS atoms within r_c of PMMA', 'SDS compounds within r_c of PMMA',
    'MMA atoms within r_c of PMMA', 'MMA compounds within r_c of PMMA',

    'PS-SDS contacts within r_c',
    'PS-MMA contacts within r_c',
    'PS-PMMA contacts within r_c',
    'PS-SOL contacts within r_c',
    'PS-Ion contacts within r_c',
    'Shell PMMA-SDS contacts within r_c',
    'Shell PMMA-MMA contacts within r_c',
    'PMMA-SOL contacts within r_c',
    'PMMA-Ion contacts within r_c',
]

compounds = ['SDS', 'MMA', 'PMMA', 'SOL', 'Na', 'PS' ]
number_of_compounds  = { 'SDS' : 79, 'MMA' : 100, 'PMMA' :  11, 'SOL' : 86223, 'Na' : 210, 'PS' : 1}
n_atoms_per_compound = { 'SDS' : 42, 'MMA' :  15, 'PMMA' : 152, 'SOL' : 3, 'Na' : 1, 'PS' : 1602 }

time_mda = data_mda[:, 2] / 1000.0 # from ps to ns

# colors for histogram
# viridis, tab10, Accent, Set1, Set2, Dark2, rainbow, turbo, jet
colors = cm.get_cmap('tab10', len(compounds) + 1)

n_running = 21
n_shift = int((n_running - 1) / 2)


##################################################################
### Histogram
# distribution of the number of SDS atoms/residues within r_c of PS, and a version normalized by the number of SDS atoms
# distribution of the number of MMA atoms/residues within r_c of PS, and a version normalized by the number of MMA atoms
# distribution of the number of PMMA atoms/residues within r_c of PS, and a version normalized by the number of PMMA atoms

if False:
    idx_atm = [
        headers.index('SDS atoms within r_c of PS'),
        headers.index('MMA atoms within r_c of PS'),
        headers.index('PMMA atoms within r_c of PS')
        ]
    idx_res = [
        headers.index('SDS compounds within r_c of PS'),
        headers.index('MMA compounds within r_c of PS'),
        headers.index('PMMA compounds within r_c of PS')
        ]
    mols = ['SDS', 'MMA', 'PMMA']

    y_shift = 750
    nbins = 10

    fig = plt.figure(tight_layout = True)
    fig.set_dpi(200)
    fig.set_size_inches(9, 6)
    gs = gridspec.GridSpec(3, 1)

    # Number of atoms
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Compounds within ' + rcs + ' of PS surface')
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[y_shift:, idx_atm[i]]
        ax.hist(y, bins = nbins, color = c, alpha = 0.3)
    ax.set_xlabel('Number of atoms')
    ax.legend()

    # Number of compounds
    ax = fig.add_subplot(gs[1, 0])
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_res[i]]
        ax.hist(y, bins = nbins, color = c, alpha = 0.3)
    ax.set_xlabel('Number of compounds')

    # Normalized number of atoms
    ax = fig.add_subplot(gs[2, 0])
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = np.zeros(data_mda[:, 1].shape)
        y[:] = data_mda[:, idx_atm[i]]
        for j in range(0, len(y)):
            if data_mda[j, idx_res[i]] > 0:
                y[j] = y[j] / (data_mda[j, idx_res[i]] * n_atoms_per_compound[compounds[i]] )
        ax.hist(y, bins = nbins, color = c, alpha = 0.3)
    ax.set_xlabel('Normalized number of atoms')
    #ax.legend()

    plt.show()


##################################################################
### Compounds within r_c of PS
# evolution of the number of SDS atoms/residues within r_c of PS, and a version normalized by the total number of atom of these SDS residues
# evolution of the number of MMA atoms/residues within r_c of PS, and a version normalized by the total number of atom of these MMA residues
# evolution of the number of PMMA atoms/residues within r_c of PS, and a version normalized by the total number of atom of these PMMA residues

if True:
    idx_atm = [
        headers.index('SDS atoms within r_c of PS'),
        headers.index('MMA atoms within r_c of PS'),
        headers.index('PMMA atoms within r_c of PS')
        ]
    idx_res = [
        headers.index('SDS compounds within r_c of PS'),
        headers.index('MMA compounds within r_c of PS'),
        headers.index('PMMA compounds within r_c of PS')
        ]
    mols = ['SDS', 'MMA', 'PMMA']

    fig = plt.figure(tight_layout = True)
    fig.set_dpi(200)
    fig.set_size_inches(9, 6)
    gs = gridspec.GridSpec(2, 1)

    # Number of atoms
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Number of atoms that are in contact with PS, $r_c=$' + rcs)
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_atm[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    #ax.set_xlabel('Time, ns')
    ax.set_ylabel('Number of atoms')
    ax.legend()
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Number of compounds
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title('Number of compounds with atoms that are in contact with PS, $r_c=$' + rcs)
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_res[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    ax.set_ylabel('Number of compounds')
    #plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Normalized number of atoms
    # Normalized by the total number of atoms of compounds on the surface - this shows the ratio of the compound atoms that are within r_c of PS surface to the total number of atom of those compounds
    # ax = fig.add_subplot(gs[2, 0])
    # ax.set_title('Normalized number of atoms that are in contact with PS, $r_c=$' + rcs)
    # for i in range(0, len(mols)):
    #     c = colors(compounds.index(mols[i]))
    #     y = np.zeros(data_mda[:, 1].shape)
    #     y[:] = data_mda[:, idx_atm[i]]
    #     for j in range(0, len(y)):
    #         if data_mda[j, idx_res[i]] > 0:
    #             y[j] = y[j] / (data_mda[j, idx_res[i]] * n_atoms_per_compound[compounds[i]] )
    #     ax.plot(time_mda, y, color = c, alpha = 0.3)
    #     y_running_mean = running_mean(y, n_running)
    #     ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    # ax.set_ylabel('Normalized number of atoms')
    # #ax.legend()

    ax.set_xlabel('Time, ns')

    plt.show()


##################################################################
### PMMA compounds within r_c of surface SDS/MMA compounds and PS
# evolution of the number of PMMA atoms/residues within r_c of surface SDS, and a version normalized by the number of SDS atoms
# evolution of the number of PMMA atoms/residues within r_c of surface MMA, and a version normalized by the number of MMA atoms
# evolution of the number of PMMA atoms/residues within r_c of PS, and a version normalized by the number of PS atoms

if True:
    idx_atm = [
        headers.index('PMMA atoms within r_c of surface SDS'),
        headers.index('PMMA atoms within r_c of surface MMA'),
            headers.index('PMMA atoms within r_c of PS')
        ]
    idx_res = [
        headers.index('PMMA compounds within r_c of surface SDS'),
        headers.index('PMMA compounds within r_c of surface MMA'),
        headers.index('PMMA compounds within r_c of PS')
        ]
    mols = ['SDS', 'MMA', 'PS']

    fig = plt.figure(tight_layout = True)
    fig.set_dpi(200)
    fig.set_size_inches(9, 6)
    gs = gridspec.GridSpec(2, 1)

    # Number of atoms
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Number of PMMA atoms in contact with PS and SDS, MMA on the PS surface, $r_c=$' + rcs)
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_atm[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    #ax.set_xlabel('Time, ns')
    ax.set_ylabel('Number of atoms')
    ax.legend()
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Number of PMMA compounds
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title('Number of PMMA compounds in contact with PS and SDS, MMA on the PS surface, $r_c=$' + rcs)
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_res[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    ax.set_ylabel('Number of compounds')
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Normalized number of atoms
    # Normalized by the total number of atoms in PMMA compounds that are within r_c of surface SDS/MMA or PS depending on the case - this shows the ratio of PMMA atoms that are in contact with other surface compounds or PS surface
    # ax = fig.add_subplot(gs[2, 0])
    # ax.set_title('... normalized by the total number of atoms in PMMA compounds that are within ' + rcs+ ' of surface SDS/MMA or PS, $r_c=$' + rcs)
    # # SDS/MMA
    # for i in range(0, len(mols)):
    #     c = colors(compounds.index(mols[i]))
    #     y = np.zeros(data_mda[:, 1].shape)
    #     y[:] = data_mda[:, idx_atm[i]]
    #     for j in range(0, len(y)):
    #         if data_mda[j, idx_res[i]] > 0:
    #             y[j] = y[j] / (data_mda[j, idx_res[i]] * n_atoms_per_compound['PMMA'] )
    #     ax.plot(time_mda, y, color = c, alpha = 0.3)
    #     y_running_mean = running_mean(y, n_running)
    #     ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    # ax.set_ylabel('Normalized number of atoms')

    ax.set_xlabel('Time, ns')

    plt.show()


##################################################################
### surface SDS/MMA compounds within r_c of PMMA compounds
if True:
    idx_atm = [
        headers.index('SDS atoms within r_c of PMMA'),
        headers.index('MMA atoms within r_c of PMMA')
        ]
    idx_res = [
        headers.index('SDS compounds within r_c of PMMA'),
        headers.index('MMA compounds within r_c of PMMA')
        ]
    mols = ['SDS', 'MMA']

    fig = plt.figure(tight_layout = True)
    fig.set_dpi(200)
    fig.set_size_inches(9, 6)
    gs = gridspec.GridSpec(2, 1)

    # Number of atoms
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Number of SDS, MMA atoms in contact with PMMA, $r_c=$' + rcs)
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_atm[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    #ax.set_xlabel('Time, ns')
    ax.set_ylabel('Number of atoms')
    ax.legend()
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Number of compounds
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title('Number of SDS, MMA compounds in contact with PMMA, $r_c=$' + rcs)
    for i in range(0, len(mols)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx_res[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    ax.set_ylabel('Number of compounds')
    # plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Normalized number of atoms
    # Normalized by the total number of atoms of compounds on the surface - this shows the ratio of the compound atoms that are within r_c of PS surface to the total number of atom of those compounds
    # ax = fig.add_subplot(gs[2, 0])
    # ax.set_title('Number of atoms in surface MMA/SDS within ' + rcs + ' of PMMA normalized by the number of atoms in these compounds, $r_c=$' + rcs)
    # for i in range(0, len(mols)):
    #     c = colors(compounds.index(mols[i]))
    #     y = np.zeros(data_mda[:, 1].shape)
    #     y = data_mda[:, idx_atm[i]]
    #     for j in range(0, len(y)):
    #         if data_mda[j, idx_res[i]] > 0:
    #             y[j] = y[j] / (data_mda[j, idx_res[i]] * n_atoms_per_compound[mols[i]] )
    #     ax.plot(time_mda, y, color = c, alpha = 0.3)
    #     y_running_mean = running_mean(y, n_running)
    #     ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    # ax.set_ylabel('Normalized number of atoms')
    # #ax.legend()

    ax.set_xlabel('Time, ns')

    plt.show()


##################################################################
### Number of PS contacts

if True:
    idx = [
        headers.index('PS-SDS contacts within r_c'),
        headers.index('PS-MMA contacts within r_c'),
        headers.index('PS-PMMA contacts within r_c'),
        headers.index('PS-SOL contacts within r_c'),
        headers.index('PS-Ion contacts within r_c')
        ]
    mols = ['SDS', 'MMA', 'PMMA', 'SOL', 'Na']

    fig = plt.figure(tight_layout = True)
    fig.set_dpi(200)
    fig.set_size_inches(9, 6)
    gs = gridspec.GridSpec(2, 1)

    # Number of atoms
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Number of contacts with PS, $r_c=$' + rcs)
    for i in range(0, len(idx)):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    #ax.set_xlabel('Time, ns')
    ax.set_ylabel('Number of contacts')
    ax.legend()
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Normalized number of contacts
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title('Ratio of the number of contacts with PS, $r_c=$' + rcs)
    total_ncontacts = np.zeros(data_mda[:, 1].shape)
    #total_ncontacts_1 = np.zeros(data_mda[:, 1].shape)
    for i in range(0, len(idx)): total_ncontacts += data_mda[:, idx[i]]
    for i in range(0, len(idx)):
        c = colors(compounds.index(mols[i]))
        y = np.zeros(data_mda[:, 1].shape)
        y[:] = data_mda[:, idx[i]]
        for j in range(0, len(y)): y[j] = y[j] / total_ncontacts[j]
        #total_ncontacts_1 += y
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c)
    #ax.plot(time_mda, total_ncontacts_1, color = 'r')
    ax.set_ylim([-0.05, 1.05])
    ax.set_ylabel('Number of contacts ratio')
    #ax.legend()

    ax.set_xlabel('Time, ns')

    plt.show()

##################################################################
### Number of PMMA contacts

if True:
    idx = [
        headers.index('Shell PMMA-SDS contacts within r_c'),
        headers.index('Shell PMMA-MMA contacts within r_c'),
        headers.index('PS-PMMA contacts within r_c'),
        headers.index('PMMA-SOL contacts within r_c'),
        headers.index('PMMA-Ion contacts within r_c'),
        ]
    include_solvent_ions = True
    len_idx = len(idx)
    if include_solvent_ions == False: len_idx = 3
    mols = ['SDS', 'MMA', 'PS', 'SOL', 'Na']

    fig = plt.figure(tight_layout = True)
    fig.set_dpi(200)
    fig.set_size_inches(9, 6)
    gs = gridspec.GridSpec(3, 1)

    # Number of atoms
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Number of contacts with shell PMMA, $r_c=$' + rcs)
    for i in range(0, len_idx):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    # ax.set_xlim([45, 305])
    #ax.set_xlabel('Time, ns')
    ax.set_ylabel('Number of contacts')
    ax.legend()
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Number of atoms
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title('Number of contacts with shell PMMA, $r_c=$' + rcs)
    for i in range(0, len(idx)-2):
        c = colors(compounds.index(mols[i]))
        y = data_mda[:, idx[i]]
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    # ax.set_xlim([45, 305])
    #ax.set_xlabel('Time, ns')
    ax.set_ylabel('Number of contacts')
    ax.legend()
    plt.setp(ax.get_xticklabels(), visible=False)   # hide x-ticks labels

    # Normalized number of contacts
    ax = fig.add_subplot(gs[2, 0])
    ax.set_title('Ratio of the number of contacts with shell PMMA, $r_c=$' + rcs)
    total_ncontacts = np.zeros(data_mda[:, 1].shape)
    #total_ncontacts_1 = np.zeros(data_mda[:, 1].shape)
    for i in range(0, len(idx)): total_ncontacts += data_mda[:, idx[i]]
    for i in range(0, len_idx):
        c = colors(compounds.index(mols[i]))
        y = np.zeros(data_mda[:, 1].shape)
        y[:] = data_mda[:, idx[i]]
        for j in range(0, len(y)):
            if total_ncontacts[j] > 0: y[j] = y[j] / total_ncontacts[j]
        #total_ncontacts_1 += y
        ax.plot(time_mda, y, color = c, alpha = 0.3)
        y_running_mean = running_mean(y, n_running)
        ax.plot(time_mda[n_shift:-n_shift], y_running_mean, color = c, label = mols[i])
    #ax.plot(time_mda, total_ncontacts_1, color = 'r')
    # ax.set_ylim([-0.05, 1.05])
    # ax.set_xlim([45, 305])
    ax.set_ylabel('Number of contacts ratio')
    ax.legend()

    ax.set_xlabel('Time, ns')

    plt.show()
