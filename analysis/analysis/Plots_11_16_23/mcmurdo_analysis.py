# MCMurdo Analysis 
# UCSC 2023
# Liam Branch, Michelle Pichardo, Robert Johnson

# Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Globals
PARAMS = {'A':4317.5, 'B':0.0545, 'n':1.5691}
# Others options
# PARAMS.update({'A':50000, 'n':2.55 })
PARAMS.update({'n':1.6298})

c = 29.9792 # cm/ns
# pmt locations
xPMT4=9.5*np.cos(np.radians(110))*2.54
yPMT4=9.5*np.sin(np.radians(110))*2.54
xPMT1=8.*np.cos(np.radians(-45))*2.54 # x2PMT1=8.*np.cos(np.radians(-53.72))*2.54 For test
yPMT1=8.*np.sin(np.radians(-45))*2.54


# load data
prefix = 'Data4ToFCalibration_'
thresh_str = ['1212_NL4151_','1616_NL4150_']
default_pattern_recog = 'RK.txt' # |xL0| < 4 cm
# without cut should increase the numbers and/or quality of the track with T1, T4 trigger
recog_without_position_cut = 'RKToF.txt' # no cut on hit position in the first layer
df = pd.read_csv(prefix+thresh_str[0]+default_pattern_recog, sep=' ')
df2 = pd.read_csv(prefix+thresh_str[1]+default_pattern_recog, sep=' ')
# print(df)

# Threshold of 12 distance from pmt
dfd1 = np.sqrt(np.power(df['x1']-xPMT1,2)+np.power(df['y1']-yPMT1,2))
dfd4 = np.sqrt(np.power(df['x4']-xPMT4,2)+np.power(df['y4']-yPMT4,2))
# Threshold of 16
df2d1 = np.sqrt(np.power(df['x1']-xPMT1,2)+np.power(df['y1']-yPMT1,2))
df2d4 = np.sqrt(np.power(df['x4']-xPMT4,2)+np.power(df['y4']-yPMT4,2))


def tof_correction(t1, t4, d1, d4, tof):
    return tof + PARAMS['A']/(1+PARAMS['B']*np.sqrt(t1)) \
         - PARAMS['A']/(1+PARAMS['B']*np.sqrt(t4))\
         + PARAMS['n']*(d4-d1) / c
         
def tof_correction_v2(d1,d4,tof):
    return tof + PARAMS['n']*(d4-d1) / c

# df_new_tof = tof_correction(df['T1'], df['T4'], dfd1, dfd4, df['ToF'])
# df2_new_tof = tof_correction(df2['T1'], df2['T4'], df2d1, df2d4, df2['ToF'])
df_new_tof = tof_correction_v2(dfd1, dfd4, df['ToF'])
df2_new_tof = tof_correction_v2(df2d1, df2d4, df2['ToF'])


fig0, ax0 = plt.subplots()
binning = np.linspace(-10,10,200)
# Threshold of 12
ax0.hist(df['ToF'], bins=binning, histtype='step', facecolor='none', edgecolor='k', label='12 ToF')
ax0.hist(df_new_tof, bins=binning, alpha=0.5, label='12 Corrected')
# Threshold of 16
ax0.hist(df2['ToF'], bins=binning, histtype='step', facecolor='none', edgecolor='k', label='16 ToF')
ax0.hist(df2_new_tof, bins=binning, alpha=0.5, label='16 Corrected')
ax0.legend()
ax0.set_ylabel('Count')
ax0.set_xlabel('ToF (ns)')
plt.show()