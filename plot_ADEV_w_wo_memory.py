#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:56:09 2022

A versatile water vapor generation module for vapor isotope
calibration and liquid isotope measurements

by Steen-Larsen and Zannoni

Atmospheric Measurement Techniques

Code for reproducing Figure 3

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
import allantools
import datetime

from import_picarro_raw_data import import_picarro_raw_data

from scipy.stats import linregress

#%%Style
from matplotlib import style as mstyle
mstyle.use('article.mplstyle')

colors = {'H2O'         :'black',
          'd17O'        :'dodgerblue',
          'd18O'        :'forestgreen',
          'dD'          :'chocolate',
          'd-excess'    :'black',
          'O17Excess'   :'black',}


#%% Load Picarro data
paths = [
        '../DATA/HKDS2156/2022/10/10/',
        '../DATA/HKDS2156/2022/10/11/',
        '../DATA/HKDS2156/2022/10/12/',
        '../DATA/HKDS2156/2022/10/13/',
        '../DATA/HKDS2156/2022/10/14/'
         ]# Directory for daily data
# Define time region
#start_date_str      = '08-08-2022 08:30:00'
#stop_date_str       = '08-08-2022 09:50:00'


# Settings
sample_oi = "BER17k"
#sample_oi = "BER17k_withmemory"
tau_oi = 3600 # seconds

analysis_type = 'freq'

calibrate_data      = True
d17O_cal_params     = [0.958390535, -0.461053701]
d18O_cal_params     = [0.96478044, -0.375228501]
dD_cal_params       = [0.937748326, -6.009660321]


df_memory = pd.read_csv('ADEV_BER17k_withmemory_CALIBRATED_R1.csv', index_col = 0)
show_memory_effect = False



#%% Load metadata for analysis
metadata_file_name = 'Timings_Picarro.xlsx'
metadata = pd.read_excel(io=metadata_file_name, sheet_name='AllanDev')
print(metadata.head(5))  # print first 5 rows of the dataframe


#%% List files and import
for file_path in paths:
    if 'Picarro_data' in locals():
        new_Picarro_data = import_picarro_raw_data(file_path, '17O')
        Picarro_data = pd.concat([Picarro_data, new_Picarro_data])
        del(new_Picarro_data)
    else:
        Picarro_data = import_picarro_raw_data(file_path, '17O')

#%% Calibrate data R1
if calibrate_data:
    print('Calibration....')
    Picarro_data['Delta_17_16'] = Picarro_data['Delta_17_16']*d17O_cal_params[0] + d17O_cal_params[1]
    Picarro_data['Delta_18_16'] = Picarro_data['Delta_18_16']*d18O_cal_params[0] + d18O_cal_params[1]
    Picarro_data['Delta_D_H'] = Picarro_data['Delta_D_H']*dD_cal_params[0] + dD_cal_params[1]
    df_memory = pd.read_csv('ADEV_BER17k_withmemory_CALIBRATED_R1.csv', index_col = 0)
    
#%% Run analysis per line
seed = 1234

for sample in metadata.index:
    if metadata['ID'][sample] == sample_oi:
        start_winOI = metadata['Avg_start'][sample]
        stop_winOI  = metadata['Avg_stop'][sample]
        Picarro_date_mask = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)
        Picarro_data_subset = Picarro_data[Picarro_date_mask].copy()
        #curr_rate = 1/(1e-9*np.mean(np.diff(Picarro_data_subset.index))/np.timedelta64(1, 'ns'))
        curr_rate = 1/(np.mean(np.diff(Picarro_data_subset.index[:3600]))/np.timedelta64(1, 's'))
        #curr_rate = 1
        # Add 17O excess variable
        Picarro_data_subset['17OExcess_calc'] = 1e6*(np.log(Picarro_data_subset['Delta_17_16']/1000 + 1) - 0.528*np.log(Picarro_data_subset['Delta_18_16']/1000 + 1))
        Picarro_data_subset['dexcess'] = Picarro_data_subset['Delta_D_H'] - 8*Picarro_data_subset['Delta_18_16']
        #taus_level          = np.arange(1, 1e4, 10)#[1, 10, 20, 60, 90, 180, 270, 360, 600, 1800, 3600, 5000]
        taus_level          = np.logspace(0, 6, 150)
        #taus_level          = np.array([1, 10, 60, 270, 1800])
        taus_level_noise    = np.logspace(0, 6, 150)       
        # Add taus of interest
        taus_level = np.concatenate((taus_level, np.array([tau_oi])))
        # Computation for H2O ------------
        (tau_out_H2O, adev_H2O, adeverr_H2O, n_H2O) = allantools.adev(Picarro_data_subset['H2O'].to_numpy(), 
                                                                               rate = curr_rate,
                                                                               data_type = analysis_type,
                                                                               taus = taus_level)  
        np.random.seed(seed)
        white_noise_H2O = np.random.randn(4*len(Picarro_data_subset['H2O']))*Picarro_data_subset['H2O'][-600:].std()+ Picarro_data_subset['H2O'][-600:].mean()
        (tau_out_H2O_noise, adev_H2O_noise, adeverr_H2O_noise, n_H2O_noise) = allantools.adev(white_noise_H2O, 
                                                                                                       rate = curr_rate,
                                                                                                       data_type = "freq",
                                                                                                       taus = taus_level_noise)
        
        # Computation for d17O ------------
        (tau_out_17_16, adev_17_16, adeverr_17_16, n_17_16) = allantools.adev(Picarro_data_subset['Delta_17_16'].to_numpy(), 
                                                                               rate = curr_rate,
                                                                               data_type = analysis_type,
                                                                               taus = taus_level)  
        np.random.seed(seed)
        white_noise_17_16 = np.random.randn(4*len(Picarro_data_subset['Delta_17_16']))*Picarro_data_subset['Delta_17_16'][-600:].std()+ Picarro_data_subset['Delta_17_16'][-600:].mean()
        (tau_out_17_16_noise, adev_17_16_noise, adeverr_17_16_noise, n_17_16_noise) = allantools.adev(white_noise_17_16, 
                                                                                                       rate = curr_rate,
                                                                                                       data_type = "freq",
                                                                                                       taus = taus_level_noise)
        
        # Computation for d18O ------------    
        (tau_out_18_16, adev_18_16, adeverr_18_16, n_18_16) = allantools.adev(Picarro_data_subset['Delta_18_16'].to_numpy(), 
                                                                               rate = curr_rate,
                                                                               data_type = analysis_type,
                                                                               taus = taus_level)
        np.random.seed(seed)
        white_noise_18_16 = np.random.randn(4*len(Picarro_data_subset['Delta_18_16']))*Picarro_data_subset['Delta_18_16'][-600:].std()+ Picarro_data_subset['Delta_18_16'][-600:].mean()
        (tau_out_18_16_noise, adev_18_16_noise, adeverr_18_16_noise, n_18_16_noise) = allantools.adev(white_noise_18_16, 
                                                                                                       rate = curr_rate,
                                                                                                       data_type = analysis_type,
                                                                                                       taus = taus_level_noise)
        # Computation for dD ------------ 
        (tau_out_D_H, adev_D_H, adeverr_D_H, n_D_H) = allantools.adev(Picarro_data_subset['Delta_D_H'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = analysis_type,
                                                                       taus = taus_level)
        
        np.random.seed(seed)
        white_noise_D_H = np.random.randn(4*len(Picarro_data_subset['Delta_D_H']))*Picarro_data_subset['Delta_D_H'][-600:].std()+ Picarro_data_subset['Delta_D_H'][-600:].mean()
        (tau_out_D_H_noise, adev_D_H_noise, adeverr_D_H_noise, n_D_H_noise) = allantools.adev(white_noise_D_H, 
                                                                                                       rate = curr_rate,
                                                                                                       data_type = analysis_type,
                                                                                                       taus = taus_level_noise)
        
        
        # Computation for 17X_calc calc ------------ 
        (tau_out_17X_calc, adev_17X_calc, adeverr_17X_calc, n_17X_calc) = allantools.adev(Picarro_data_subset['17OExcess_calc'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = analysis_type,
                                                                       taus = taus_level)
        np.random.seed(seed)
        white_noise_17X_calc = np.random.randn(4*len(Picarro_data_subset['17OExcess_calc']))*Picarro_data_subset['17OExcess_calc'][-600:].std()+ Picarro_data_subset['17OExcess_calc'][-600:].mean()
        (tau_out_17X_calc_noise, adev_17X_calc_noise, adeverr_17X_calc_noise, n_17X_calc_noise) = allantools.adev(white_noise_17X_calc, 
                                                                                                       rate = curr_rate,
                                                                                                       data_type = analysis_type,
                                                                                                       taus = taus_level_noise)
        
        # Computation for 17X from Picarro ------------ 
        (tau_out_17X, adev_17X, adeverr_17X, n_17X) = allantools.adev(1000*Picarro_data_subset['Excess_17'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = analysis_type,
                                                                       taus = taus_level)
        np.random.seed(seed)
        white_noise_17X = np.random.randn(4*len(Picarro_data_subset['Excess_17']))*(1000*Picarro_data_subset['Excess_17'][-600:]).std()+ (1000*Picarro_data_subset['Excess_17'][-600:]).mean()
        (tau_out_17X_noise, adev_17X_noise, adeverr_17X_noise, n_17X_noise) = allantools.adev(white_noise_17X, 
                                                                                               rate = curr_rate,
                                                                                               data_type = analysis_type,
                                                                                               taus = taus_level_noise)
        # Computation for d-excess        ------------ 
        (tau_out_dX, adev_dX, adeverr_dX, n_dX) = allantools.adev(Picarro_data_subset['dexcess'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = analysis_type,
                                                                       taus = taus_level)
        np.random.seed(seed)
        white_noise_dX = np.random.randn(4*len(Picarro_data_subset['dexcess']))*(Picarro_data_subset['dexcess'][-600:]).std()+ (Picarro_data_subset['dexcess'][-600:]).mean()
        (tau_out_dX_noise, adev_dX_noise, adeverr_dX_noise, n_dX_noise) = allantools.adev(white_noise_dX, 
                                                                                               rate = curr_rate,
                                                                                               data_type = analysis_type,
                                                                                               taus = taus_level_noise)
        # Computation for cavity temperature        ------------ 
        (tau_out_CT, adev_CT, adeverr_CT, n_CT) = allantools.adev(Picarro_data_subset['CavityTemp'].to_numpy(), 
                                                                  rate = curr_rate,
                                                                  data_type = analysis_type,
                                                                  taus = taus_level)
        # Computation for cavity pressure        ------------ 
        (tau_out_CP, adev_CP, adeverr_CP, n_CP) = allantools.adev(Picarro_data_subset['CavityPressure'].to_numpy(), 
                                                                  rate = curr_rate,
                                                                  data_type = "freq",
                                                                  taus = taus_level)
    
#%% 
length = np.timedelta64(stop_winOI-start_winOI)/np.timedelta64(1, 's')
buff_title = "Sample: %s (@%d ± %d ppmv) - Duration: %d s" % (sample_oi, Picarro_data_subset['H2O'].mean(), Picarro_data_subset['H2O'].std(), length)

#%% Load injection timins
df_injections = pd.read_csv('SP_BER_inj_time.csv')
df_injections['Start'] = pd.to_datetime(df_injections['Start'], dayfirst = True)
df_injections['Stop'] = pd.to_datetime(df_injections['Stop'], dayfirst = True)
df_injections.drop(index = 13, inplace = True)


#%% Subplots
window_AVG = 270

fig, ax = plt.subplots(nrows = 3, figsize = [4,5], dpi = 300,
                       gridspec_kw={'height_ratios': [0.5, 1, 1]})
# Water
ax[0].plot(Picarro_data_subset.index,
           Picarro_data_subset['H2O'].rolling(window_AVG).mean(),
           color = colors['H2O'])
ax[0].plot(Picarro_data.index,
           Picarro_data['H2O'].rolling(window_AVG).mean(),
           color = colors['H2O'],
           alpha = 0.2)           
ax[0].set_ylabel('H$_2$O (ppmv)')
ax[0].set_ylim([15000, 20000])
ax[0].set_xlim([pd.to_datetime('10/10/2022 16:00', dayfirst=True), pd.to_datetime('14/10/2022 11:30', dayfirst=True), ])
ax[0].grid('on', which='both')
ax[0].set_xticklabels([])
ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%d-%H'))
ax[0].set_xlabel('Time (day-hour)')
# Panel label
ax[0].text(0.95, 0.2, '(a)', horizontalalignment='center',
     verticalalignment='center', transform=ax[0].transAxes, weight = 1000, fontsize=14)

top0, bottom0 = ax[0].get_ylim()
booly = True
for inj in range(df_injections['Start'].count()):
    ax[0].text(df_injections['Start'][inj], 18300, df_injections['T'][inj])
    if booly:
            if df_injections['T'][inj] == 1:
                curr_col = [1,1,1]#'dodgerblue'
            elif df_injections['T'][inj] == 2:
                curr_col = [.75, 0.75,.75]#'palegreen'
            elif df_injections['T'][inj] == 3:
                curr_col = [.5,.5,.5]#'plum'
            elif df_injections['T'][inj] == 4:
                curr_col = curr_col = [.25,.25,.25]#'salmon'
            #curr_col = 'dodgerblue'            
            ax[0].fill_betweenx([bottom0, top0],
                                df_injections['Start'][inj], df_injections['Stop'][inj],
                                color = curr_col, alpha=0.25)


# Isotopes, permil
ax[1].loglog(tau_out_17_16, adev_17_16, label = '$\delta^{17}$O', 
             color = colors['d17O'])
if show_memory_effect:
    ax[1].loglog(df_memory.index, df_memory['d17O'], 
                 color = colors['d17O'], ls = '--')

ax[1].loglog(tau_out_18_16, adev_18_16, label = '$\delta^{18}$O', 
             color = colors['d18O'])
if show_memory_effect:
    ax[1].loglog(df_memory.index, df_memory['d18O'], 
                 color = colors['d18O'], ls = '--')

ax[1].loglog(tau_out_D_H, adev_D_H, label = '$\delta$D', 
             color = colors['dD'])
if show_memory_effect:
    ax[1].loglog(df_memory.index, df_memory['dD'],
                 color = colors['dD'], ls = '--')

ax[1].loglog(tau_out_dX, adev_dX, label = 'd-excess', 
             color = colors['d-excess'])
if show_memory_effect:
    ax[1].loglog(df_memory.index, df_memory['dX'], 
                 color = colors['d-excess'], ls = '--')

#R1 ADD NOISE
ax[1].loglog(tau_out_18_16_noise[tau_out_18_16_noise < 1e3], adev_18_16_noise[tau_out_18_16_noise < 1e3], label = 'White noise', 
             color = 'r', lw = 1, ls = ':')

ax[1].set_ylim([1e-3, 2])
ax[1].set_ylabel('$\delta$ (‰)')
ax[1].legend(fontsize = 5, loc = 'lower left')
ax[1].set_xticklabels([])
# Panel label
ax[1].text(0.95, 0.85, '(b)', horizontalalignment='center',
     verticalalignment='center', transform=ax[1].transAxes, weight = 1000, fontsize=14)

# Isotopes, permeg
ax[2].loglog(tau_out_17X_calc, adev_17X_calc, label = '$\Delta^{17}$O', 
             color = colors['O17Excess'])
if show_memory_effect:
    ax[2].loglog(df_memory.index, df_memory['X17O'], 
                 color = colors['O17Excess'], ls = '--')

#R1 ADD NOISE
ax[2].loglog(tau_out_17X_calc_noise[tau_out_17X_calc_noise < 1e3], adev_17X_calc_noise[tau_out_17X_calc_noise < 1e3], label = 'White noise', 
             color = 'r', lw = 1, ls = ':')


ax[2].set_ylabel('$\Delta^{17}$O (per meg)')

ax[2].set_xlabel('Time (s)')

ax[2].legend(fontsize = 5, loc = 'lower left')

# Panel label
ax[2].text(0.95, 0.85, '(c)', horizontalalignment='center',
     verticalalignment='center', transform=ax[2].transAxes, weight = 1000, fontsize=14)

for axis in ax:
    axis.grid('on', which = 'both')

plt.tight_layout(h_pad = 0)

fig.savefig('asd.pdf')