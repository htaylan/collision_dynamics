import pandas as pd
import itertools
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import os
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import time as time
from multiprocessing import Process, Pool, Manager

nproc = 16 # The number of cores to use
num_of_files = 175 #The number of trajectory files
times = [100,200,500,1000,1500,1750] #The life-time criteria for formed complex
par_distances = [4,5,7.5,10,12.5] #The distance criteria for formed complex 
cut_off = 20 #If distance larger than this value, stop calculating
manager = Manager()
matrix_values = manager.list()

input_directory ='cis_dcd_imp0' #The directory of dcd files
output_directory='cis_dist_imp0' #The directory to write distance files
psf_file = 'cis_dbb_ca1.psf' #The psf file of the system

starttime = time.time()


def calc_com_br(i):
    #Provide the path and the filename
    dcd_file = os.path.join(os.getcwd(), input_directory,'dbbca_test'+str(i)+'.dcd')
    Universe = mda.Universe(psf_file,dcd_file)
    time_dist = []

    for frame in Universe.trajectory:
            #Atom selection for distance calculation
            Ca  = Universe.select_atoms("name CA")
            BrC = Universe.select_atoms("name BR or name BR1").center_of_mass()
            dist = distances.distance_array(BrC,
                                            Ca.positions,
                                            box=Universe.dimensions)
            if dist >= cut_off:
                break
            else:
                lines = "{:.2f} {:.5f}".format(Universe.trajectory.time, float(dist))
                time_dist.append(lines)

        #Writes out distance files, provide the path
    dist_file  = os.path.join(os.getcwd(), output_directory, 'dist_com_'+str(i)+'.dat')
    with open(dist_file,'w') as f:
        for t_dist in time_dist:
            f.writelines("{}\n".format(t_dist))

combination = list(itertools.product(range(len(par_distances)),range(len(times))))

def analyze_ltime_dist(k):
    i, j = combination[k]
    dist_to_br = par_distances[i]
    formation_time = times[j]
    formed_complexes = []
    for i_files in range(1,num_of_files+1):
        try:
            fname = os.path.join(os.getcwd(), output_directory, 'dist_com_'+str(i_files)+'.dat')
            df =  pd.read_csv(fname,header=None,sep='\s+',usecols=[1],names='D')
            liste = df.index[df['D']<=dist_to_br].tolist()
            count = 0
            for i_list in range(len(liste)-1):
                if liste[i_list]+1 == liste[i_list+1]:
                    count += 0.1
                    if count >= formation_time:
                         formed_complexes.append(fname)
                         break
                else:
                    count = 0
        except:
            pass


    matrix_values.append([dist_to_br,formation_time,len(formed_complexes)])

    #Writes out the .dcd names which met the criteria
    out_file_name = 'D'+str(dist_to_br) + '_' + 'T'+str(formation_time) + '_'+str(len(formed_complexes))+'_imp0.dat'
    with open(out_file_name,'w') as f:
        for i_newfile in range(len(formed_complexes)):
            f.write(str(i+1)+' '+formed_complexes[i_newfile]+'\n')



with Pool(nproc) as p:
    p.map(calc_com_br, range(1,num_of_files+1))

with Pool(nproc) as p:
    p.map(analyze_ltime_dist, range(len(combination)))



matrix_values = sorted(matrix_values, key=lambda x: (int(x[0]), x[1]))
matrix_values = [x[2] for x in matrix_values]


#Plot the heat map
matrix_values = np.asarray(matrix_values).reshape(len(par_distances),len(times))

fig, ax = plt.subplots()
im = ax.imshow(matrix_values)


ax.set_xticks(np.arange(len(times)))
ax.set_yticks(np.arange(len(par_distances)))

ax.set_xticklabels(times)
ax.set_yticklabels(par_distances)

ax.set_xlabel('Time (ps)')
ax.set_ylabel('Distance (A)')

plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")


for i in range(len(par_distances)):
    for j in range(len(times)):
        text = ax.text(j, i, matrix_values[i, j],
                       ha="center", va="center", color="w")

ax.set_title("Number of Complex Formation")
fig.tight_layout()
#plt.show()
#Change the title of the plot
plt.savefig('Heat_map_longer_time_com2.png')
print('Analysis took {} seconds'.format(time.time() - starttime))
