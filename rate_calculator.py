import pandas as pd
import numpy as np
import argparse

"""
Example usage as follows:

    python3 rate_calculator.py -f filename.dat -w Stratified

"""

# Constants used in the rate calculations
dbb_mass = 211.88216 # Unit: au
ca_mass = 40.078 # Unit: au
kb = 1.38064852 * 10**-23 # Unit: kgm^2/s^2K 
T = 10 # Unit: K 
m_to_cm = 100 # Convert meter to cm
au_to_kg = 10 ** -27 # Convert au to kg 

# Get the filename and other arguments 
parser = argparse.ArgumentParser('Calculate Rate')
parser.add_argument('--filename', '-f', help='Filename')
parser.add_argument('--wmethod', '-w', choices=['weight','Stratified'], default='weight')
args = parser.parse_args()

# Calculate the reduced mass of the reactants in kg
def rmass_calc(m1,m2):
    m1 = m1 * 1.6605 * au_to_kg
    m2 = m2 * 1.6605 * au_to_kg
    rmass = m1 * m2 / (m1+m2)
    return rmass

# Calculate the weight for all b values
def weight_calc(b,bmax):
    weight = 2*b / bmax
    return weight

# Strafied weight
def sweight_calc(b_max,b_min,bmax):
    sweight = (b_max**2 - b_min**2) / bmax**2
    return sweight

# Calculate the probability of complex formation at the given b value
def preac_calc(weight,prob):
    preac = np.sum(weight * prob)
    return preac

# Calculate the rate of the reaction
def rate_calc():
    rate = np.sqrt(8*kb*T/(np.pi*rmass)) * m_to_cm *  np.pi * bmax_rate**2 * preac
    return rate

##### Rate calculations #####

# Load the b and Prob values as a dataframe
df = pd.read_csv(args.filename,header=None,sep='\s+')

# Give name to the columns
df.columns=['b','Prob']

# Get the bmax value in angstorm
bmax = df.shape[0]

#Convert bmax value to cm 
bmax_rate = bmax * 10 ** -8

# Calculate the weight of each b based on different weight methods. 
if args.wmethod == 'weight':
    df['Weight'] = weight_calc(df['b'],bmax)


if args.wmethod == 'Stratified':

    if df['b'][0] - 1 >= 0 :
        df['Weight'] = sweight_calc(df['b'],df['b']-1,bmax)

    elif df['b'][0] -1 < 0:
        df['Weight'] = sweight_calc(df['b'],df['b']-1,bmax)
        df['Weight'][0] = 0


# Calculate the rate
preac = preac_calc(df.Weight,df.Prob)
rmass = rmass_calc(dbb_mass,ca_mass)
rate = rate_calc()
print(rate)
