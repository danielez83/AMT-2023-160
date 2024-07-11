#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:56:09 2022

A versatile water vapor generation module for vapor isotope
calibration and liquid isotope measurements

by Steen-Larsen and Zannoni

Atmospheric Measurement Techniques

Code for reproducing Figure 8

Update the "paths" list with the actual path of PICARRO data

@author: daniele.zannoni@uib.no
"""

#%% Add path to Python
import os
import sys
#cur_dir = os.path.dirname(__file__)
#sys.path.append(os.path.join(cur_dir, "/functions"))

#%% Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import ticker as mticker
import allantools
import datetime

from import_picarro_raw_data import import_picarro_raw_data

from scipy.stats import linregress

from scipy import signal

from numpy.fft import rfft, rfftfreq

#%%Style
from matplotlib import style as mstyle
mstyle.use('article.mplstyle')

colors = {'H2O'         :'black',
          'd17O'        :'dodgerblue',
          'd18O'        :'forestgreen',
          'dD'          :'chocolate',
          'd-excess'    :'black',
          'O17Excess'   :'black',}

style_2092 = "--"
style_2156 = "-"

#%% Load Picarro data
paths_2156 = [
            #'../Picarro Data/2022/10/10/',
            '../DATA/HKDS2156/2022/05/09/',
            '../DATA/HKDS2156/2022/05/10/',
            '../DATA/HKDS2156/2022/05/11/',
            ]# Directory for daily data
paths_2092 = [
            #'../Picarro Data/2022/10/10/',
            '../DATA/HKDS2092/2022/05/09/',
            '../DATA/HKDS2092/2022/05/10/',
            '../DATA/HKDS2092/2022/05/11/',
            ]# Directory for daily data
# Define time region
#start_date_str      = '08-08-2022 08:30:00'
#stop_date_str       = '08-08-2022 09:50:00'


# Settings
#sample_oi = "BER14k"
#sample_oi = "AMB"
sample_oi = ["BER14k",      
             "AMB"]         # 5 hours

calibrate_humidity = True


#%% Load metadata for analysis
metadata_file_name = 'Timings_Picarro.xlsx'
metadata = pd.read_excel(io=metadata_file_name, sheet_name='Picarros_comparison')
print(metadata.head(5))  # print first 5 rows of the dataframe


#%% List files and import for 2092
for file_path in paths_2092:
    if 'Picarro_data_2092' in locals():
        new_Picarro_data_2092 = import_picarro_raw_data(file_path, '17O')
        Picarro_data_2092 = pd.concat([Picarro_data_2092, new_Picarro_data_2092])
        del(new_Picarro_data_2092)
    else:
        Picarro_data_2092 = import_picarro_raw_data(file_path, '17O')
#%% List files and import for 2156
for file_path in paths_2156:
    if 'Picarro_data_2156' in locals():
        new_Picarro_data_2156 = import_picarro_raw_data(file_path, '17O')
        Picarro_data_2156 = pd.concat([Picarro_data_2156, new_Picarro_data_2156])
        del(new_Picarro_data_2156)
    else:
        Picarro_data_2156 = import_picarro_raw_data(file_path, '17O')
        
#%% Run analysis per line
seed = 1234

for sample in metadata.index:
    if metadata['ID'][sample] in str(sample_oi):
        print("Loading --> " + metadata['ID'][sample])
        start_winOI = metadata['Avg_start'][sample]
        stop_winOI  = metadata['Avg_stop'][sample]
        
        Picarro_2092_date_mask = (Picarro_data_2092.index > start_winOI) & (Picarro_data_2092.index < stop_winOI)
        Picarro_data_2092_subset = Picarro_data_2092[Picarro_2092_date_mask].copy()
        Picarro_data_2092_subset.interpolate(inplace = True)
        
        Picarro_2156_date_mask = (Picarro_data_2156.index > start_winOI) & (Picarro_data_2156.index < stop_winOI)
        Picarro_data_2156_subset = Picarro_data_2156[Picarro_2156_date_mask].copy()
        Picarro_data_2156_subset.interpolate(inplace = True)
        
        # Resample timeseries
        Picarro_data_2092_subset_resampled = Picarro_data_2092_subset.resample('1S').mean()
        Picarro_data_2092_subset_resampled.interpolate(inplace = True)
        Picarro_data_2156_subset_resampled = Picarro_data_2156_subset.resample('1S').mean()
        Picarro_data_2156_subset_resampled.interpolate(inplace = True)
        
        # Calibrate humidity
        if calibrate_humidity:
            Picarro_data_2092_subset['H2O'] = Picarro_data_2092_subset['H2O'] * 0.733863 + 33.805885 # Use WIFVOS calibration based on VAISALA HMP113 data (22/04/10)
            Picarro_data_2092_subset_resampled['H2O'] = Picarro_data_2092_subset_resampled['H2O'] * 0.733863 + 33.805885 # Use WIFVOS calibration based on VAISALA HMP113 data (22/04/10)
            # Calibrate 2156 humidity with 2092 data
            #slope, intercept, r, pval, stderr = linregress(Picarro_data_2156_subset_resampled['H2O'], Picarro_data_2092_subset_resampled['H2O'])
            #Picarro_data_2156_subset['H2O'] = Picarro_data_2156_subset['H2O'] * slope + intercept
            #if sample_oi == 'AMB':
                #Picarro_data_2156_subset['H2O'] = Picarro_data_2156_subset['H2O'] * 1.0017101928465322 + 39.58942970428143
                #Picarro_data_2156_subset_resampled['H2O'] = Picarro_data_2156_subset_resampled['H2O'] * 1.0017101928465322 + 39.58942970428143
            #else:
                #Picarro_data_2156_subset['H2O'] = Picarro_data_2156_subset['H2O'] * 0.9437834540624788 + 879.6438668803567
                #Picarro_data_2156_subset_resampled['H2O'] = Picarro_data_2156_subset_resampled['H2O'] * 0.9437834540624788 + 879.6438668803567
        # Save Picarro data in dataframe
        strbuff = "Picarro_data_2092_subset"+"_"+metadata['ID'][sample]+"=Picarro_data_2092_subset.copy()"
        exec(strbuff)
        strbuff = "Picarro_data_2092_subset_resampled"+"_"+metadata['ID'][sample]+"=Picarro_data_2092_subset_resampled.copy()"
        exec(strbuff)
        strbuff = "Picarro_data_2156_subset"+"_"+metadata['ID'][sample]+"=Picarro_data_2156_subset.copy()"
        exec(strbuff)
        strbuff = "Picarro_data_2156_subset_resampled"+"_"+metadata['ID'][sample]+"=Picarro_data_2156_subset_resampled.copy()"
        exec(strbuff)

#%% Show humidity stats for BER14k and AMB
# fig, ax = plt.subplots()
# ax.plot(np.arange(0, len(Picarro_data_2092_subset_resampled_AMB)), 
#         Picarro_data_2092_subset_resampled_AMB['H2O'], 
#         label = '2092 Ambient')
# ax.plot(np.arange(0, len(Picarro_data_2092_subset_resampled_BER14k)), 
#         Picarro_data_2092_subset_resampled_BER14k['H2O'], 
#         label = '2092 Cal. Dev.')
# ax.plot(np.arange(0, len(Picarro_data_2156_subset_resampled_AMB)), 
#         Picarro_data_2156_subset_resampled_AMB['H2O'], 
#         label = '2156 Ambient')
# ax.plot(np.arange(0, len(Picarro_data_2156_subset_resampled_BER14k)), 
#         Picarro_data_2156_subset_resampled_BER14k['H2O'], 
#         label = '2156 Cal. Dev.')
# ax.legend()
print("Length of BER14k dataframe is: %d" % len(Picarro_data_2092_subset_resampled_BER14k))
print("Length of AMB dataframe is: %d" % len(Picarro_data_2092_subset_resampled_AMB))
print("Mean +/- 1 sigma for BER14k with HKDS 2156 is: %d +/- %d" % (Picarro_data_2156_subset_BER14k['H2O'].mean(), Picarro_data_2156_subset_BER14k['H2O'].std()))
print("Mean +/- 1 sigma for AMB with HKDS 2156 is: %d +/- %d" % (Picarro_data_2156_subset_AMB['H2O'].mean(), Picarro_data_2156_subset_AMB['H2O'].std()))
mdl = linregress(np.arange(0, len(Picarro_data_2156_subset_resampled_BER14k)), Picarro_data_2156_subset_resampled_BER14k['H2O'])
print("Trend for BER14k with HKDS 2156 is: %d ppm/h" % (mdl.slope*3600))
mdl = linregress(np.arange(0, len(Picarro_data_2156_subset_resampled_AMB)), Picarro_data_2156_subset_resampled_AMB['H2O'])
print("Trend for AMB with HKDS 2156 is: %d ppm/h" % (mdl.slope*3600))
print("The two timeseries are comparable in term of variability (different average value!)")

#%% Print correlations
print("Correlations between 2092 and 2156")
print("H2O " + str(np.corrcoef(Picarro_data_2092_subset_resampled_BER14k['H2O'],
                              Picarro_data_2156_subset_resampled_BER14k['H2O'])[0][1]))
print("d17O " + str(np.corrcoef(Picarro_data_2092_subset_resampled_BER14k['Delta_17_16'], 
                                Picarro_data_2156_subset_resampled_BER14k['Delta_17_16'])[0][1]))
print("d18O " + str(np.corrcoef(Picarro_data_2092_subset_resampled_BER14k['Delta_18_16'], 
                                Picarro_data_2156_subset_resampled_BER14k['Delta_18_16'])[0][1]))
print("dD " + str(np.corrcoef(Picarro_data_2092_subset_resampled_BER14k['Delta_D_H'], 
                              Picarro_data_2156_subset_resampled_BER14k['Delta_D_H'])[0][1]))
#%% Calculate FFT of H2O and isotopes
# For H2O the FFT of original timeseries and resampled timeseries are nearly identical
# For d18O the FFT of original timeseries and resampled timeseries are nearly identical
# For d17O the FFT of original timeseries and resampled timeseries are nearly identical, 
# slight differences for high freq > 1^-1 Hz
# For dDO the FFT of original timeseries and resampled timeseries are nearly identical, 
# slight differences for high freq > 1^-1 Hz
# Can't do the same analysis for Delta of BER14k vs AMB, H2O and Delta are highly correlated in natural water vapor
#%% FFT for H2O, ambient and cal dev.
# Perform FFT analysis of H2O signal 2092 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_BER14k['H2O'].count()
# H2O --------------------------------        
f_H2O_2092_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_BER14k = rfft(Picarro_data_2092_subset_BER14k['H2O'].values)
P_H2O_2092_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_BER14k)**2)
# Perform FFT analysis of H2O signal 2092 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_AMB['H2O'].count()
# H2O --------------------------------        
f_H2O_2092_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_AMB = rfft(Picarro_data_2092_subset_AMB['H2O'].values)
P_H2O_2092_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_AMB)**2)
# Perform FFT analysis of H2O signal 2156 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_BER14k['H2O'].count()
# H2O --------------------------------        
f_H2O_2156_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_BER14k = rfft(Picarro_data_2156_subset_BER14k['H2O'].values)
P_H2O_2156_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_BER14k)**2)
# Perform FFT analysis of H2O signal 2156 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_AMB['H2O'].count()
# H2O --------------------------------        
f_H2O_2156_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_AMB = rfft(Picarro_data_2156_subset_AMB['H2O'].values)
P_H2O_2156_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_AMB)**2)

#%% FFT for Delta_18_16, ambient and cal dev.
# Perform FFT analysis of Delta_18_16 signal 2092 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_BER14k['Delta_18_16'].count()
# Delta_18_16 --------------------------------        
f_Delta_18_16_2092_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_BER14k = rfft(Picarro_data_2092_subset_BER14k['Delta_18_16'].values)
P_Delta_18_16_2092_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_BER14k)**2)
# Perform FFT analysis of Delta_18_16 signal 2092 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_AMB['Delta_18_16'].count()
# Delta_18_16 --------------------------------        
f_Delta_18_16_2092_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_AMB = rfft(Picarro_data_2092_subset_AMB['Delta_18_16'].values)
P_Delta_18_16_2092_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_AMB)**2)
# Perform FFT analysis of Delta_18_16 signal 2156 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_BER14k['Delta_18_16'].count()
# Delta_18_16 --------------------------------        
f_Delta_18_16_2156_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_BER14k = rfft(Picarro_data_2156_subset_BER14k['Delta_18_16'].values)
P_Delta_18_16_2156_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_BER14k)**2)
# Perform FFT analysis of Delta_18_16 signal 2156 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_AMB['Delta_18_16'].count()
# Delta_18_16 --------------------------------        
f_Delta_18_16_2156_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_AMB = rfft(Picarro_data_2156_subset_AMB['Delta_18_16'].values)
P_Delta_18_16_2156_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_AMB)**2)

#%% FFT for Delta_17_16, ambient and cal dev.
# Perform FFT analysis of Delta_17_16 signal 2092 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_BER14k['Delta_17_16'].count()
# Delta_17_16 --------------------------------        
f_Delta_17_16_2092_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_BER14k = rfft(Picarro_data_2092_subset_BER14k['Delta_17_16'].values)
P_Delta_17_16_2092_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_BER14k)**2)
# Perform FFT analysis of Delta_17_16 signal 2092 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_AMB['Delta_17_16'].count()
# Delta_17_16 --------------------------------        
f_Delta_17_16_2092_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_AMB = rfft(Picarro_data_2092_subset_AMB['Delta_17_16'].values)
P_Delta_17_16_2092_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_AMB)**2)
# Perform FFT analysis of Delta_17_16 signal 2156 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_BER14k['Delta_17_16'].count()
# Delta_17_16 --------------------------------        
f_Delta_17_16_2156_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_BER14k = rfft(Picarro_data_2156_subset_BER14k['Delta_17_16'].values)
P_Delta_17_16_2156_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_BER14k)**2)
# Perform FFT analysis of Delta_17_16 signal 2156 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_AMB['Delta_17_16'].count()
# Delta_17_16 --------------------------------        
f_Delta_17_16_2156_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_AMB = rfft(Picarro_data_2156_subset_AMB['Delta_17_16'].values)
P_Delta_17_16_2156_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_AMB)**2)

#%% FFT for Delta_D_H, ambient and cal dev.
# Perform FFT analysis of Delta_D_H signal 2092 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_BER14k['Delta_D_H'].count()
# Delta_D_H --------------------------------        
f_Delta_D_H_2092_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_BER14k = rfft(Picarro_data_2092_subset_BER14k['Delta_D_H'].values)
P_Delta_D_H_2092_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_BER14k)**2)
# Perform FFT analysis of Delta_D_H signal 2092 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_AMB['Delta_D_H'].count()
# Delta_D_H --------------------------------        
f_Delta_D_H_2092_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_AMB = rfft(Picarro_data_2092_subset_AMB['Delta_D_H'].values)
P_Delta_D_H_2092_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_AMB)**2)
# Perform FFT analysis of Delta_D_H signal 2156 Calibration device
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_BER14k['Delta_D_H'].count()
# Delta_D_H --------------------------------        
f_Delta_D_H_2156_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_BER14k = rfft(Picarro_data_2156_subset_BER14k['Delta_D_H'].values)
P_Delta_D_H_2156_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_BER14k)**2)
# Perform FFT analysis of Delta_D_H signal 2156 Ambient
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_AMB.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_AMB['Delta_D_H'].count()
# Delta_D_H --------------------------------        
f_Delta_D_H_2156_AMB = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_AMB = rfft(Picarro_data_2156_subset_AMB['Delta_D_H'].values)
P_Delta_D_H_2156_AMB = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_AMB)**2)

#%% FFT for Cavity temperature
# 2092
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_BER14k['CavityTemp'].count()
# H2O --------------------------------        
f_CT_2092_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_BER14k = rfft(Picarro_data_2092_subset_BER14k['CavityTemp'].values)
P_CT_2092_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_BER14k)**2)
# 2156
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_BER14k['CavityTemp'].count()
# CT --------------------------------        
f_CT_2156_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_BER14k = rfft(Picarro_data_2156_subset_BER14k['CavityTemp'].values)
P_CT_2156_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_BER14k)**2)

#%% FFT for Cavity pressure
# 2092
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2092_subset_BER14k['CavityPressure'].count()
# --------------------------------        
f_CP_2092_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2092_BER14k = rfft(Picarro_data_2092_subset_BER14k['CavityPressure'].values)
P_CP_2092_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2092_BER14k)**2)
# 2156
# Sampling frequency
fs = 1/(np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
N = Picarro_data_2156_subset_BER14k['CavityPressure'].count()
# CP --------------------------------        
f_CP_2156_BER14k = rfftfreq(N, 1/fs)
#Autpower, scale autopower accounting for 2-sides of the spectrum and converto to PSD
y_hat_2156_BER14k = rfft(Picarro_data_2156_subset_BER14k['CavityPressure'].values)
P_CP_2156_BER14k = 2 * (1/(fs*N)) * (np.abs(y_hat_2156_BER14k)**2)



#%% Allan deviation
taus_level = np.logspace(0, 5, 150)
taus_level_noise    = np.logspace(0, 6, 150) 
curr_rate_2092 = (np.mean(np.diff(Picarro_data_2092_subset_BER14k.index))/np.timedelta64(1, 's'))
# Computation for d17O ------------
(tau_out_17_16_2092, adev_17_16_2092, adeverr_17_16_2092, n_17_16_2092) = allantools.adev(Picarro_data_2092_subset['Delta_17_16'].to_numpy(), 
                                                                                          rate = curr_rate_2092,
                                                                                          data_type = "freq",
                                                                                          taus = taus_level)  
# Computation for d18O ------------    
(tau_out_18_16_2092, adev_18_16_2092, adeverr_18_16_2092, n_18_16_2092) = allantools.adev(Picarro_data_2092_subset['Delta_18_16'].to_numpy(), 
                                                                                          rate = curr_rate_2092,
                                                                                          data_type = "freq",
                                                                                          taus = taus_level)

# Computation for dD ------------ 
(tau_out_D_H_2092, adev_D_H_2092, adeverr_D_H_2092, n_D_H_2092) = allantools.adev(Picarro_data_2092_subset['Delta_D_H'].to_numpy(), 
                                                                                  rate = curr_rate_2092,
                                                                                  data_type = "freq",
                                                                                  taus = taus_level)

# Computation for Cavity Temperature ------------ 
(tau_out_CT_2092, adev_CT_2092, adeverr_CT_2092, n_CT_2092) = allantools.adev(Picarro_data_2092_subset['CavityTemp'].to_numpy(), 
                                                                              rate = curr_rate_2092,
                                                                              data_type = "freq",
                                                                              taus = taus_level)

# Computation for Cavity Pressure ------------ 
(tau_out_CP_2092, adev_CP_2092, adeverr_CP_2092, n_CP_2092) = allantools.adev(Picarro_data_2092_subset['CavityPressure'].to_numpy(), 
                                                                              rate = curr_rate_2092,
                                                                              data_type = "freq",
                                                                              taus = taus_level)


curr_rate_2156 = (np.mean(np.diff(Picarro_data_2156_subset_BER14k.index))/np.timedelta64(1, 's'))
# Computation for d17O ------------
(tau_out_17_16_2156, adev_17_16_2156, adeverr_17_16_2156, n_17_16_2156) = allantools.adev(Picarro_data_2156_subset['Delta_17_16'].to_numpy(), 
                                                                                          rate = curr_rate_2156,
                                                                                          data_type = "freq",
                                                                                          taus = taus_level)  
# Computation for d18O ------------    
(tau_out_18_16_2156, adev_18_16_2156, adeverr_18_16_2156, n_18_16_2156) = allantools.adev(Picarro_data_2156_subset['Delta_18_16'].to_numpy(), 
                                                                                          rate = curr_rate_2156,
                                                                                          data_type = "freq",
                                                                                          taus = taus_level)
np.random.seed(seed)
white_noise_18_16 = np.random.randn(4*len(Picarro_data_2156_subset['Delta_18_16']))*Picarro_data_2156_subset['Delta_18_16'][-600:].std()+ Picarro_data_2156_subset['Delta_18_16'][-600:].mean()
(tau_out_18_16_noise, adev_18_16_noise, adeverr_18_16_noise, n_18_16_noise) = allantools.adev(white_noise_18_16, 
                                                                                               rate = curr_rate_2092,
                                                                                               data_type = 'freq',
                                                                                               taus = taus_level_noise)

# Computation for dD ------------ 
(tau_out_D_H_2156, adev_D_H_2156, adeverr_D_H_2156, n_D_H_2156) = allantools.adev(Picarro_data_2156_subset['Delta_D_H'].to_numpy(), 
                                                                                  rate = curr_rate_2156,
                                                                                  data_type = "freq",
                                                                                  taus = taus_level)

# Computation for Cavity Temperature ------------ 
(tau_out_CT_2156, adev_CT_2156, adeverr_CT_2156, n_CT_2156) = allantools.adev(Picarro_data_2156_subset['CavityTemp'].to_numpy(), 
                                                                              rate = curr_rate_2156,
                                                                              data_type = "freq",
                                                                              taus = taus_level)

# Computation for Cavity Pressure ------------ 
(tau_out_CP_2156, adev_CP_2156, adeverr_CP_2156, n_CP_2156) = allantools.adev(Picarro_data_2156_subset['CavityPressure'].to_numpy(), 
                                                                              rate = curr_rate_2156,
                                                                              data_type = "freq",
                                                                              taus = taus_level)

#%% Plot timeseries
#smooth timeseries for better visualization
# 300 is the maximum in Allan deviation of Cavity Temperature
avg_time = 300
# Moving average filter -------------------------------------------------------
avg_time = 8 * 60
Picarro_data_2092_subset_smooth = Picarro_data_2092_subset.rolling(avg_time).mean()
Picarro_data_2156_subset_smooth = Picarro_data_2156_subset.rolling(avg_time).mean()
# Butterworth filter ----------------------------------------------------------
#Picarro_data_2092_subset_smooth = Picarro_data_2092_subset.copy()
#Picarro_data_2156_subset_smooth = Picarro_data_2156_subset.copy()
#omega = 1/(10*60)
#b,a = signal.butter(1, omega, 'lowpass')
#for column in Picarro_data_2092_subset.columns:
#    buff_data = signal.filtfilt(b,a, Picarro_data_2092_subset[column])
#    Picarro_data_2092_subset_smooth[column] = buff_data
#    buff_data = signal.filtfilt(b,a, Picarro_data_2156_subset[column])
#    Picarro_data_2156_subset_smooth[column] = buff_data

# Resampling and synch
Picarro_data_2092_subset_smooth_resampled = Picarro_data_2092_subset_smooth.resample('1S').mean()
Picarro_data_2156_subset_smooth_resampled = Picarro_data_2156_subset_smooth.resample('1S').mean()
Picarro_data_2092_subset_smooth_resampled.interpolate(inplace = True)
Picarro_data_2156_subset_smooth_resampled.interpolate(inplace = True)

# Data must be synched, estimated delay is 1288 seconds
# 2092 delyaed respect to 2156
Picarro_data_2092_subset_smooth_resampled = Picarro_data_2092_subset_smooth_resampled.iloc[1000:-1288, :]
Picarro_data_2156_subset_smooth_resampled = Picarro_data_2156_subset_smooth_resampled.iloc[1288:-1000, :]

# Normalization of useful data
d18O_2092_norm  = (Picarro_data_2092_subset_smooth_resampled['Delta_18_16'] - Picarro_data_2092_subset_smooth_resampled['Delta_18_16'].mean())/Picarro_data_2092_subset_smooth_resampled['Delta_18_16'].std()
d18O_2156_norm  = (Picarro_data_2156_subset_smooth_resampled['Delta_18_16'] - Picarro_data_2156_subset_smooth_resampled['Delta_18_16'].mean())/Picarro_data_2156_subset_smooth_resampled['Delta_18_16'].std()
d18O_2092       = Picarro_data_2092_subset_smooth_resampled['Delta_18_16']
d18O_2156       = Picarro_data_2156_subset_smooth_resampled['Delta_18_16']
dD_2092_norm    = (Picarro_data_2092_subset_smooth_resampled['Delta_D_H'] - Picarro_data_2092_subset_smooth_resampled['Delta_D_H'].mean())/Picarro_data_2092_subset_smooth_resampled['Delta_D_H'].std()
dD_2156_norm    = (Picarro_data_2156_subset_smooth_resampled['Delta_D_H'] - Picarro_data_2156_subset_smooth_resampled['Delta_D_H'].mean())/Picarro_data_2156_subset_smooth_resampled['Delta_D_H'].std()
d17O_2092_norm  = (Picarro_data_2092_subset_smooth_resampled['Delta_17_16'] - Picarro_data_2092_subset_smooth_resampled['Delta_17_16'].mean())/Picarro_data_2092_subset_smooth_resampled['Delta_17_16'].std()
d17O_2156_norm  = (Picarro_data_2156_subset_smooth_resampled['Delta_17_16'] - Picarro_data_2156_subset_smooth_resampled['Delta_17_16'].mean())/Picarro_data_2156_subset_smooth_resampled['Delta_17_16'].std()
CT_2092_norm    = (Picarro_data_2092_subset_smooth_resampled['CavityTemp'] - Picarro_data_2092_subset_smooth_resampled['CavityTemp'].mean())/Picarro_data_2092_subset_smooth_resampled['CavityTemp'].std()
CT_2156_norm    = (Picarro_data_2156_subset_smooth_resampled['CavityTemp'] - Picarro_data_2156_subset_smooth_resampled['CavityTemp'].mean())/Picarro_data_2156_subset_smooth_resampled['CavityTemp'].std()
CT_2092         = Picarro_data_2092_subset_smooth_resampled['CavityTemp']
CT_2156         = Picarro_data_2156_subset_smooth_resampled['CavityTemp']


#%% Create single figure
#import matplotlib.gridspec as gridspec

#fig = plt.figure(tight_layout=True, figsize=(6,3))
#gs = gridspec.GridSpec(1, 2)

#ax_ADEV = fig.add_subplot(gs[:, 0])
#ax_d18O = fig.add_subplot(gs[0, 1])
#ax_CT = fig.add_subplot(gs[1, 1])

fig, (ax_ADEV, ax_CT) = plt.subplots(ncols = 2, figsize=(5,2.5))

#---
ax_ADEV.loglog(tau_out_17_16_2092, adev_17_16_2092, color = colors['d17O'], linestyle = style_2092)#, label = "HKDS2092 - $\delta^{17}$O")
ax_ADEV.loglog(tau_out_18_16_2092, adev_18_16_2092, color = colors['d18O'], linestyle = style_2092)#, label = "HKDS2092 - $\delta^{18}$O")
ax_ADEV.loglog(tau_out_D_H_2092, adev_D_H_2092, color = colors['dD'], linestyle = style_2092)#, label = "HKDS2092 - $\delta$D")
ax_ADEV.loglog(tau_out_17_16_2156, adev_17_16_2156, color = colors['d17O'], linestyle = style_2156, label = "$\delta^{17}$O")#label = "HKDS2156 - $\delta^{17}$O")
ax_ADEV.loglog(tau_out_18_16_2156, adev_18_16_2156, color = colors['d18O'], linestyle = style_2156, label = "$\delta^{18}$O")#label = "HKDS2156 - $\delta^{18}$O")
ax_ADEV.loglog(tau_out_D_H_2156, adev_D_H_2156, color = colors['dD'], linestyle = style_2156, label = "$\delta$D")#label = "HKDS2156 - $\delta$D")

ax_ADEV.loglog(tau_out_18_16_noise[tau_out_18_16_noise < 3e2], adev_18_16_noise[tau_out_18_16_noise < 3e2],
               label = 'White noise', color = 'r', lw = 1, ls = ':')

ax_ADEV.set_ylim([3e-3, 1e0])
ax_ADEV.set_ylabel('$\delta$ (‰)')
ax_ADEV.set_xlabel('Time (s)')
ax_ADEV.legend(loc = 'lower left', fontsize = "xx-small")
ax_ADEV.xaxis.set_major_locator(mticker.LogLocator(numticks=999))
ax_ADEV.xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
ax_ADEV.grid(which = 'minor')

# Panel label
ax_ADEV.text(.90, .90, '(a)', horizontalalignment='center',
             verticalalignment='center', transform=ax_ADEV.transAxes, weight = 1000, fontsize=14)

#---
limits = [-3.5, 3.5]
ax_CT.scatter(CT_2092_norm, d18O_2092_norm, marker = 'o', color = 'k', facecolor = 'w', alpha = .1, s = 4)
model_2092 = linregress(CT_2092_norm, d18O_2092_norm)
strbuff = "R$_{2092}$ = %.2f" % model_2092.rvalue
ax_CT.text(-3.3, 2.8, strbuff, color = 'k')


ax_CT.scatter(CT_2156_norm, d18O_2156_norm, marker = 'x', color = 'r',alpha = .1, s = 4)
model_2156 = linregress(CT_2156_norm, d18O_2156_norm)
strbuff = "R$_{2156}$ = %.2f" % model_2156.rvalue
ax_CT.text(-3.3, 2.3, strbuff, color = 'r')

ax_CT.set_ylabel('$\delta^{18}$O$_{Z}$')
ax_CT.set_xlabel('Cavity T$_{Z}$')
ax_CT.set_xlim(limits)
ax_CT.set_ylim(limits)
ax_CT.grid(which = 'minor')

# Panel label
ax_CT.text(.90, .10, '(b)', horizontalalignment='center',
           verticalalignment='center', transform=ax_CT.transAxes, weight = 1000, fontsize=14)

plt.tight_layout()
#fig.savefig('asd.pdf')

#%%
