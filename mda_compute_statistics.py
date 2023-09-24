# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse
import time

import MDAnalysis as mda
from MDAnalysis.analysis import contacts

def contacts_within_cutoff(ts, group_a, group_b, r_cutoff = 2.5):
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions)
    # determine which distances <= r_cutoff
    n_contacts = contacts.contact_matrix(dist, r_cutoff).sum()
    return n_contacts

#process inputs
parser = argparse.ArgumentParser()
parser.add_argument("-b", help="begin time (ps)", default=000000, type=float)
parser.add_argument("-e", help="end time (ps)", default=300000, type=float)
parser.add_argument("-dt" ,help="time step in trajectory (ps)", default=200, type=float)
parser.add_argument("-dta", help="time step for analysis (ps)", default=200, type=float)
args = parser.parse_args()

begTime = args.b
endTime = args.e
dt = args.dt
dta = args.dta

# radius of cutoff [angstroms]
r_cutoff = 2.5 # angstroms
r_cutoff_shell = 10 # angstroms
# string representation of r_cutoff to use in selection commands
rc = str(r_cutoff)
rshell = str(r_cutoff_shell)

folder = "./20211104-part3_pollyMMA-300ns/"
filename = folder + "INIT_POLYMMA_prepared_em_nvt_npt_md_5B6centered"
outputFilenamePrefix = 'stats_mda_frames_' + str(int(begTime/dt)) + '-' + str(int(endTime/dt)) + '_' + rc + 'A'

u = mda.Universe(filename+".tpr", filename+".xtc", in_memory=False)

#########
# system:
# 5B6 = Polystyrene (PS)
# AZM = SDS
# WVD = MMA
# 6CB = PMMA
#########

# Atom selection language in MDAnalysis is used to perform selections
# see: https://userguide.mdanalysis.org/stable/selections.html

# atoms in various compounds
ps = u.select_atoms("resname 5B6")
pmma = u.select_atoms("resname 6CB")
sds = u.select_atoms("resname AZM")
mma = u.select_atoms("resname WVD")
sol = u.select_atoms("resname SOL")
ions = u.select_atoms("resname NA")
non_ps = u.select_atoms("resname 6CB or resname AZM or resname WVD or resname SOL or resname NA")

pmma_atoms = pmma.atoms

# atoms in PMMA residues that surround PS
# NSL atoms:                   (n.t a in "6CB") within 3A of ("5B6")
# NSL residues: n.t sg having ((n.t a in "6CB") within 3A of ("5B6"))
pmma_around_ps = pmma.select_atoms('around ' + rc + ' global group ps', updating=True, ps=ps)
pmma_shell_ps = pmma.select_atoms('around ' + rshell + ' global group ps', updating=True, ps=ps)
# atoms in SDS residues that surround PS
# NSL atoms:                   (n.t a in "AZM") within 3A of ("5B6")
# NSL residues: n.t sg having ((n.t a in "AZM") within 3A of ("5B6"))
sds_around_ps = sds.select_atoms('around ' + rc + ' global group ps', updating=True, ps=ps)
sds_shell_ps = sds.select_atoms('around ' + rshell + ' global group ps', updating=True, ps=ps)
# atoms in MMA residues that surround PS
# NSL atoms:                   (n.t a in "WVD") within 3A of ("5B6")
# NSL residues: n.t sg having ((n.t a in "WVD") within 3A of ("5B6"))
mma_around_ps = mma.select_atoms('around ' + rc + ' global group ps', updating=True, ps=ps)
mma_shell_ps = mma.select_atoms('around ' + rshell + ' global group ps', updating=True, ps=ps)
# atoms in NA residues that surround PS
# NSL atoms:                   (n.t a in "NA") within 3A of ("5B6")
# NSL residues: n.t sg having ((n.t a in "NA") within 3A of ("5B6"))
ions_around_ps = ions.select_atoms('around ' + rc + ' global group ps', updating=True, ps=ps)
# atoms in SOL residues that surround PS
# NSL atoms:                   (n.t a in "SOL") within 3A of ("5B6")
# NSL residues: n.t sg having ((n.t a in "SOL") within 3A of ("5B6"))
sol_around_ps = sol.select_atoms('around ' + rc + ' global group ps', updating=True, ps=ps)

nframes = int((endTime - begTime) / dta) + 1
nFrame = len(u.trajectory)

df_column_headers = [
    'frame', 'time (ps)',

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

df = pd.DataFrame(index = np.arange(0, nframes), columns = df_column_headers)

#print(' | '.join(df.columns.values.tolist()))

iframe = -1
tfirst = time.time()

for ts in u.trajectory[int(begTime/dt):]:
    current_time = u.trajectory.time
    if (u.trajectory.time >= begTime) and (u.trajectory.time <= endTime) and (u.trajectory.time % dta == 0 ):
        tframe = time.time()
        iframe += 1

        # SDS around PS
        #all_atoms_of_sds_around_ps = sds_around_ps.residues.atoms
        all_atoms_of_sds_shell_ps = sds_shell_ps.residues.atoms
        # MMA around PS
        #all_atoms_of_mma_around_ps = mma_around_ps.residues.atoms
        all_atoms_of_mma_shell_ps = mma_shell_ps.residues.atoms
        # PMMA around PS
        #all_atoms_of_pmma_around_ps = pmma_around_ps.residues.atoms
        all_atoms_of_pmma_shell_ps = pmma_shell_ps.residues.atoms

        # atoms in PS that are exposed to solvent, ions, SDS, MMA, PMMA
        # NSL atoms: (n.t a in "5B6") within 2.5A of ((n.t a in "AZM") or (n.t a in "WVD") or (n.t a in "6CB") or (n.t a in "SOL") or (n.t a in "NA"))
        # ps_accessible = ps.select_atoms('around 2.5 global group non_ps', updating=True, non_ps = non_ps)
        # ps_accessible2mma = ps.select_atoms('around 2.5 global group mma', updating=True, mma=mma)
        # ps_accessible2sds = ps.select_atoms('around 2.5 global group sds', updating=True, sds=sds)
        # ps_accessible2pmma = ps.select_atoms('around 2.5 global group pmma', updating=True, pmma=pmma)
        # ps_accessible2ions = ps.select_atoms('around 2.5 global group ions', updating=True, ions=ions)
        # ps_accessible2sol = ps.select_atoms('around 2.5 global group sol', updating=True, sol=sol)
        # t = [ps_accessible.n_atoms, ps_accessible2mma.n_atoms, ps_accessible2sds.n_atoms, ps_accessible2pmma.n_atoms, ps_accessible2ions.n_atoms, ps_accessible2sol.n_atoms, ps_accessible2mma.n_atoms + ps_accessible2sds.n_atoms + ps_accessible2pmma.n_atoms + ps_accessible2ions.n_atoms + ps_accessible2sol.n_atoms]
        # print(' | '.join(str(x).rjust(4) for x in t))

        # number of contacts between PS and various groups
        #ncont_ps_nonps = contacts_within_cutoff(ts, ps, non_ps) # this is equal to the sum of the contacs per compound
        ncont_ps_sds   = contacts_within_cutoff(ts, ps, sds_around_ps, r_cutoff)
        ncont_ps_mma   = contacts_within_cutoff(ts, ps, mma_around_ps, r_cutoff)
        ncont_ps_pmma  = contacts_within_cutoff(ts, ps, pmma_around_ps, r_cutoff)
        ncont_ps_sol   = contacts_within_cutoff(ts, ps, sol_around_ps, r_cutoff)
        ncont_ps_ions  = contacts_within_cutoff(ts, ps, ions_around_ps, r_cutoff)

        # atoms in PMMA that are around SDS in shell that surrounds PS
        # NSL atoms:                   (n.t a in "6CB") within 3A of (n.t sg having ((n.t a in "AZM") within 3A of ("5B6")))
        # NSL residues: n.t sg having ((n.t a in "6CB") within 3A of (n.t sg having ((n.t a in "AZM") within 3A of ("5B6"))))
        pmma_around_sds = pmma.select_atoms('around ' + rc + ' global group sds_shell', updating=True, sds_shell = all_atoms_of_sds_shell_ps)
        all_atoms_of_pmma_around_sds = pmma_around_sds.residues.atoms

        # atoms in PMMA that are around MMA shell that surrounds PS
        # NSL atoms:                   (n.t a in "6CB") within 3A of (n.t sg having ((n.t a in "WVD") within 3A of ("5B6")))
        # NSL residues: n.t sg having ((n.t a in "6CB") within 3A of (n.t sg having ((n.t a in "WVD") within 3A of ("5B6"))))
        pmma_around_mma = pmma.select_atoms('around ' + rc + ' global group mma_shell', updating=True, mma_shell = all_atoms_of_mma_shell_ps)
        all_atoms_of_pmma_around_mma = pmma_around_mma.residues.atoms

        # atoms in surface SDS that are around PMMA that are around SDS that surround PS
        # NSL comment: atoms in SDS surrounding PS intersect with atoms in SDS within 3A of PMMA
        # NSL atoms:                     (n.t a in (n.t sg having ((n.t a in "AZM") within 3A of ("5B6")))) and ((n.t a in "AZM") within 3A of (n.t a in "6CB"))
        # NSL residues:   n.t sg having ((n.t a in (n.t sg having ((n.t a in "AZM") within 3A of ("5B6")))) and ((n.t a in "AZM") within 3A of (n.t a in "6CB")))
        sds_around_ps_pmma = all_atoms_of_sds_shell_ps.select_atoms('around ' + rc + ' global group pmma_around_sds', updating=True, pmma_around_sds = all_atoms_of_pmma_around_sds)

        # atoms in surface MMA that are around PMMA that are around SDS that surround PS
        # NSL comment: atoms in MMA surrounding PS intersect with atoms in MMA within 3A of PMMA
        # NSL atoms:                     (n.t a in (n.t sg having ((n.t a in "WVD") within 3A of ("5B6")))) and ((n.t a in "WVD") within 3A of (n.t a in "6CB"))
        # NSL residues:   n.t sg having ((n.t a in (n.t sg having ((n.t a in "WVD") within 3A of ("5B6")))) and ((n.t a in "WVD") within 3A of (n.t a in "6CB")))
        mma_around_ps_pmma = all_atoms_of_mma_shell_ps.select_atoms('around ' + rc + ' global group pmma_around_mma', updating=True, pmma_around_mma = all_atoms_of_pmma_around_mma)

        # solvent around PMMA in PS shell
        sol_around_ps_pmma = sol.select_atoms('around ' + rc + ' global group all_atoms_of_pmma_shell_ps', updating=True, all_atoms_of_pmma_shell_ps = all_atoms_of_pmma_shell_ps)
        # ions around PMMA in PS shell
        ions_around_ps_pmma = ions.select_atoms('around ' + rc + ' global group all_atoms_of_pmma_shell_ps', updating=True, all_atoms_of_pmma_shell_ps = all_atoms_of_pmma_shell_ps)

        # number of contacts between PMMA and surface SDS/MMA
        ncont_pmma_sds = contacts_within_cutoff(ts, all_atoms_of_pmma_shell_ps, all_atoms_of_sds_shell_ps, r_cutoff)
        ncont_pmma_mma = contacts_within_cutoff(ts, all_atoms_of_pmma_shell_ps, all_atoms_of_mma_shell_ps, r_cutoff)
        ncont_pmma_sol = contacts_within_cutoff(ts, all_atoms_of_pmma_shell_ps, sol_around_ps_pmma, r_cutoff)
        ncont_pmma_ions = contacts_within_cutoff(ts, all_atoms_of_pmma_shell_ps, ions_around_ps_pmma, r_cutoff)

        # fill in the data row
        row_data = [ts.frame, current_time,
                    # SDS around PS
                    sds_around_ps.n_atoms, sds_around_ps.n_residues,
                    # MMA around PS
                    mma_around_ps.n_atoms, mma_around_ps.n_residues,
                    # PMMA around PS
                    pmma_around_ps.n_atoms, pmma_around_ps.n_residues,
                    # atoms in PMMA that are around SDS that surround PS
                    pmma_around_sds.n_atoms, pmma_around_sds.n_residues,
                    # atoms in PMMA that are around MMA that surround PS
                    pmma_around_mma.n_atoms, pmma_around_mma.n_residues,
                    # atoms in SDS around PMMA that are around SDS that surround PS
                    sds_around_ps_pmma.n_atoms, sds_around_ps_pmma.n_residues,
                    # atoms in MMA around PMMA that are around SDS that surround PS
                    mma_around_ps_pmma.n_atoms, mma_around_ps_pmma.n_residues,
                    # number of contacts between PS and other compounds
                    ncont_ps_sds, ncont_ps_mma, ncont_ps_pmma, ncont_ps_sol, ncont_ps_ions,
                    # number of contacts between shell PMMA and other compounds
                    ncont_pmma_sds, ncont_pmma_mma, ncont_pmma_sol, ncont_pmma_ions
                    ]
        # add data to the data frame
        df.loc[iframe] = row_data
        print(' | '.join(str(x).rjust(4) for x in row_data))

print("Processed in {:8.4f} s".format(time.time() - tfirst))

# save the data frame to a file
df.to_csv(outputFilenamePrefix + ".csv")
