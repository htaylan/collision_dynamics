#!/usr/bin/python
#

import pandas as pd
import numpy as np
import subprocess
import random
import os
from multiprocessing import Process, Pool, Manager, cpu_count

# nproc = cpu_count() - 2  #Number of CPU to use
nproc = 1  # Number of CPU to use
imp_factor = 1  # The impact factor value
num_of_file = 1000  # The number of initial geometries to generate
des_velocity = 843  # Velocity in units of m/s
old_file_path = '/home/taylan/Desktop/dbb-ca/old_test/dbbca_crd_imp'+str(imp_factor)


def a_akma_to_m_s(velocity):
    #velocity / 100 converts m/s to A/ps, then converted to A/AKMA by 20/0.978 
    con_factor = (velocity / 100) / (20/0.978)
    return con_factor

con_factor = a_akma_to_m_s(des_velocity)

#rate (cm*3/s-1)  to lifetime (ps)
def rate_to_lifetime(rate,bmax):
    #Convert cm^3/s-1 to A^3/s-1  
    rate_A = rate  * 10 ** 24
    rate_s = rate / bmax**3
    lifetime = 1/rate_s * 10**-12
    #lifetime = (1 / (rate  * 10 ** 24 * bmax**3))
    return print(lifetime)

rate_to_lifetime(1.4*10**-9,15)

#Generate new geometries based on DBB geometry files
def norm_velo_vec(i):
    #Define file name
    filename = 'dbbca_'+str(i)+'_imp'+str(imp_factor)+'.vel'

    #Read file as dataframe, and extract x,y,z components of the velocity vector
    df = pd.read_csv(filename,skiprows=range(0,4),sep='\s+',header=None)

    #Get old velocity components, generate new y and z coordinates
    comp_x = df.iloc[10,4] ; comp_y_old = df.iloc[10,5] ; comp_z_old =  df.iloc[10,6]
    comp_z = comp_z_old + imp_factor + random.uniform(0,1)
    comp_y = np.sqrt((comp_y_old**2 - comp_z**2))
    vector = np.sqrt(comp_x**2+comp_y**2+comp_z**2)

    #Unit conversion from  m/s to A/AKMA

    norm_comp_x = (comp_x / vector * con_factor)
    norm_comp_y = -(comp_y / vector * con_factor)
    norm_comp_z = -(comp_z / vector * con_factor)

    #Coordinate file name
    cor_filename = 'dbbca_'+str(i)+'_imp'+str(imp_factor)+'.cor'

    #Replace 800 m/s velocity components with the old ones (with impact factor 10 or more, delete the space before " {:.5f}/g" on the last line
    process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_x,norm_comp_x,filename)], shell=True)
    process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_y_old,norm_comp_y,filename)], shell=True)
    process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_z_old,norm_comp_z,filename)], shell=True)
    process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_y_old,comp_y,cor_filename)], shell=True)
    process = subprocess.Popen(['sed -i "s/{:.5f}/ {:.5f}/g" {}'.format(comp_z_old,comp_z,cor_filename)], shell=True)


#Scale velocities of new initial geometries to desired velocity. Only works for y,z axis
def scale_current_velocities(range_vel):
    new_velocities = [400,450,500,550,600,650,700,750]

    for i in range(1,range_vel+1):
        for j in range(len(new_velocities)):
            filename = 'dbbca_'+str(i)+'_imp'+str(imp_factor)+'.vel'
            df = pd.read_csv(filename, skiprows=range(0,4), sep='\s+',header=None)
            comp_x = df.iloc[10,4] ; comp_y_old = df.iloc[10,5] ; comp_z_old =  df.iloc[10,6]
            scale_to_new_velo = new_velocities[j] / des_velocity

            comp_y_new = comp_y_old * scale_to_new_velo
            comp_z_new = comp_z_old * scale_to_new_velo

            new_vel_name = 'dbbca_'+str(i)+'_imp'+str(imp_factor)+'_'+str(new_velocities[j])+'.vel'

            process = subprocess.Popen(['cp {} {}'.format(filename,new_vel_name)], shell=True)
            process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_y_old,comp_y_new,new_vel_name)], shell=True)
            process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_z_old,comp_z_new,new_vel_name)], shell=True)

old_file_path = os.getcwd()

#Scale first initial geometries to desired velocity
def scale_old_to_new(i):

   old_filename = os.path.join(old_file_path,'simu_'+str(i)+'.vel')
   df = pd.read_csv(old_filename,skiprows=range(0,3),sep='\s+',header=None)

   comp_x_old = df.iloc[10,4] ; comp_y_old = df.iloc[10,5] ; comp_z_old =  df.iloc[10,6]

   vector = np.sqrt(comp_x_old**2+comp_y_old**2+comp_z_old**2)

   sfac = con_factor / vector

   comp_x, comp_y, comp_z = sfac * comp_x_old, sfac * comp_y_old, sfac * comp_z_old

   new_vector = np.sqrt(comp_x**2+comp_y**2+comp_z**2)

   process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_x_old,comp_x,old_filename)], shell=True)
   process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_y_old,comp_y,old_filename)], shell=True)
   process = subprocess.Popen(['sed -i "s/{:.5f}/{:.5f}/g" {}'.format(comp_z_old,comp_z,old_filename)], shell=True)


for i in range(1,16):
    rate_to_lifetime(3*10**-10,i)


#with Pool(nproc) as p:
#    p.map(scale_old_to_new, range(1,num_of_file+1))

#scale_old_to_new(1)

#scale_current_velocities(100)

#with Pool(nproc) as p:
#    p.map(norm_velo_vec, range(1,num_of_file+1))
