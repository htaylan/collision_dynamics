import os
import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool, Manager
#from MDAnalysis.analysis.dihedrals import Dihedral

impact_parameter = 'imp0'
input_directory ='/data/taylan/dbbca/new_velocity_800/cis/imp0' #The directory of dcd files
psf_file = 'cis_dbb_ca1.psf' #The psf file of the system
bmax = 15

def get_succesful_trajectories():
    # Get initial conditions which formed a complex
    filenames = [fname for fname in os.listdir(impact_parameter) if fname.startswith("D10")]
    return filenames[0]

def get_traj_names(filename):
    file_path = os.path.join(impact_parameter,filename)
    df = pd.read_csv(file_path,header=None, sep = ' ',usecols=[1])
    df.columns=['F']
    df['N'] = df['F'].apply(lambda x: x.split('.')[0].split('_')[-1])
    return df

def get_distances(df,i):
    df_distances = pd.read_csv(df['F'][i], header=None, usecols=[1], sep = ' ')
    df_distances.columns=['dist']
    distances = df_distances['dist'].tolist()
    return distances

def calculate_dihedral(df,i):

    # Get the trajectory filename
    traj_file = 'dbbca_'+str(df['N'][0])+'_mdcm_uncons.dcd'
    dcd_file = os.path.join(input_directory,traj_file)

    # Initialize the dihedral analysus
    Universe = mda.Universe(psf_file,dcd_file)

    # Select dihedral angle atoms
    r=Universe.select_atoms("bynum 9 4 5 10")
    dihedral_angle = r.dihedral


    # Calculate dihedral angles
    dihedral_list = []

    for frame in Universe.trajectory:
        dihedral_list.append(dihedral_angle.value())

    return dihedral_list

filename = get_succesful_trajectories()
df = get_traj_names(filename)
distances = get_distances(df,122)
dihedral_list = calculate_dihedral(df,122)

# Plot Dist vs Dihedral Plot
plt.scatter(distances,dihedral_list)
#plt.hist2d(distances,dihedral_list,18)
plt.xlabel('Distance')
plt.ylabel('Dihedral')
plt.show()
