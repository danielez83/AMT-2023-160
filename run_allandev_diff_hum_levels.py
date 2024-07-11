#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:56:09 2022

A versatile water vapor generation module for vapor isotope
calibration and liquid isotope measurements

by Steen-Larsen and Zannoni

Atmospheric Measurement Techniques

Code for reproducing Figure 5 and Figure 6

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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
        #'../DATA/HKDS2156/2022/10/10/',
        '../DATA/HKDS2156/2022/10/11/',
        '../DATA/HKDS2156/2022/10/12/',
        '../DATA/HKDS2156/2022/10/13/',
        '../DATA/HKDS2156/2022/10/14/',
        '../DATA/HKDS2156/2022/11/13/',
        '../DATA/HKDS2156/2022/11/14/',
        '../DATA/HKDS2156/2022/11/15/',
        '../DATA/HKDS2156/2022/11/18/',
        '../DATA/HKDS2156/2022/11/19/',
        '../DATA/HKDS2156/2022/11/20/',
        '../DATA/HKDS2156/2022/11/21/'
         ]# Directory for daily data
# Define time region
#start_date_str      = '08-08-2022 08:30:00'
#stop_date_str       = '08-08-2022 09:50:00'


# Settings
sample_oi = ["BER20k",
             "BER17kST_2",
             "BER16k",
             "BER13k",
             "BER11k",
             "BER11kbis",
             "BER9k",
             "BER7k",
             "BER5k",
             "BER2k4",
             "BER0k9",
             "BER0k5"]

taus_level = np.array([1, 10, 60, 120, 270, 600, 1800, 3600])  # seconds

calibrate_data      = True
d17O_cal_params     = [0.958390535, -0.461053701]
d18O_cal_params     = [0.96478044, -0.375228501]
dD_cal_params       = [0.937748326, -6.009660321]


#%% Load metadata for analysis
metadata_file_name = 'Timings_Picarro.xlsx'
metadata = pd.read_excel(io=metadata_file_name, sheet_name='AllanDevST')
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
    #df_memory = pd.read_csv('ADEV_BER17k_withmemory_CALIBRATED_R1.csv', index_col = 0)

#%% Build dataframe to store results
col_names = ['H2O', 'std', '1', '10', '60', '120', '270', '600', '1800', '3600']
df_results_H2O = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns=col_names)
df_results_d17 = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns=col_names)
df_results_d18 = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns=col_names)
df_results_dD  = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns=col_names)
df_results_dx  = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns=col_names)
df_results_17X = pd.DataFrame([[0,0,0,0,0,0,0,0,0,0]], columns=col_names)

# R1 ramp timeseries
df_ramp_H2O = pd.DataFrame({'ID':'a',
                            'H2O':0},
                           index = [0])

#%% Setup an array to store linear regression paramters
lin_regress_params = np.array([])


#%% Run analysis per line
H2O_rows_list       = []
d17O_rows_list      = []
d18O_rows_list      = []
dD_rows_list        = []
dx_rows_list        = []
D17X_calc_rows_list  = []
D17X_rows_list       = []


# https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=50505
# why oadev is better than adev. Increase of confidence
for sample in metadata.index:
    if metadata['ID'][sample] in str(sample_oi):
        start_winOI = metadata['Avg_start'][sample]
        stop_winOI  = metadata['Avg_stop'][sample]
        Picarro_date_mask = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)
        Picarro_data_subset = Picarro_data[Picarro_date_mask].copy()
        curr_rate = 1
        # Add 17O excess variable
        Picarro_data_subset['17OExcess_calc'] = 1e6*(np.log(Picarro_data_subset['Delta_17_16']/1000 + 1) - 0.528*np.log(Picarro_data_subset['Delta_18_16']/1000 + 1))
        Picarro_data_subset['dexcess'] = Picarro_data_subset['Delta_D_H'] - 8*Picarro_data_subset['Delta_18_16']
        #taus_level          = np.arange(1, 1e4, 10)#[1, 10, 20, 60, 90, 180, 270, 360, 600, 1800, 3600, 5000]
        #taus_level          = np.logspace(0, 6, 150)   
        # Computation for H2O ------------
        (tau_out_H2O, adev_H2O, adeverr_H2O, n_H2O) = allantools.oadev(Picarro_data_subset['H2O'].to_numpy(), 
                                                                               rate = curr_rate,
                                                                               data_type = "freq",
                                                                               taus = taus_level)  
        
        # Computation for d17O ------------
        (tau_out_17_16, adev_17_16, adeverr_17_16, n_17_16) = allantools.oadev(Picarro_data_subset['Delta_17_16'].to_numpy(), 
                                                                               rate = curr_rate,
                                                                               data_type = "freq",
                                                                               taus = taus_level)  
        
        # Computation for d18O ------------    
        (tau_out_18_16, adev_18_16, adeverr_18_16, n_18_16) = allantools.oadev(Picarro_data_subset['Delta_18_16'].to_numpy(), 
                                                                               rate = curr_rate,
                                                                               data_type = "freq",
                                                                               taus = taus_level)

        # Computation for dD ------------ 
        (tau_out_D_H, adev_D_H, adeverr_D_H, n_D_H) = allantools.oadev(Picarro_data_subset['Delta_D_H'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = "freq",
                                                                       taus = taus_level)
        
        
        # Computation for 17X_calc calc ------------ 
        (tau_out_17X_calc, adev_17X_calc, adeverr_17X_calc, n_17X_calc) = allantools.oadev(Picarro_data_subset['17OExcess_calc'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = "freq",
                                                                       taus = taus_level)
        
        # Computation for 17X from Picarro ------------ 
        (tau_out_17X, adev_17X, adeverr_17X, n_17X) = allantools.oadev(1000*Picarro_data_subset['Excess_17'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = "freq",
                                                                       taus = taus_level)
        # Computation for d-excess        ------------ 
        (tau_out_dX, adev_dX, adeverr_dX, n_dX) = allantools.oadev(Picarro_data_subset['dexcess'].to_numpy(), 
                                                                       rate = curr_rate,
                                                                       data_type = "freq",
                                                                       taus = taus_level)
        # Save data into dataframes H2O
        new_row = {'H2O': Picarro_data_subset['H2O'].mean(),
                   'std':Picarro_data_subset['H2O'].std(),
                   'run lenght':(Picarro_data_subset.count()[0]/curr_rate)/3600,
                   '1':     adev_H2O[0],
                   '10':    adev_H2O[1],
                   '60':    adev_H2O[2],
                   '120':   adev_H2O[3],
                   '270':   adev_H2O[4],
                   '600':   adev_H2O[5],
                   '1800':  adev_H2O[6],
                   '3600':  adev_H2O[7]}
        #df_results_H2O = df_results_H2O.append(new_row, ignore_index = True)
        H2O_rows_list.append(new_row)

        # Save data into dataframes d17O
        new_row = {'H2O': Picarro_data_subset['H2O'].mean(),
                   'std':Picarro_data_subset['H2O'].std(),
                   'run lenght':(Picarro_data_subset.count()[0]/curr_rate)/3600,
                   '1':     adev_17_16[0],
                   '10':    adev_17_16[1],
                   '60':    adev_17_16[2],
                   '120':   adev_17_16[3],
                   '270':   adev_17_16[4],
                   '600':   adev_17_16[5],
                   '1800':  adev_17_16[6],
                   '3600':  adev_17_16[7]}
        #df_results_d17 = df_results_d17.append(new_row, ignore_index = True)
        d17O_rows_list.append(new_row)
        
        # Save data into dataframes d18O
        new_row = {'H2O': Picarro_data_subset['H2O'].mean(),
                   'std':Picarro_data_subset['H2O'].std(),
                   'run lenght':(Picarro_data_subset.count()[0]/curr_rate)/3600,
                   '1':     adev_18_16[0],
                   '10':    adev_18_16[1],
                   '60':    adev_18_16[2],
                   '120':   adev_18_16[3],
                   '270':   adev_18_16[4],
                   '600':   adev_18_16[5],
                   '1800':  adev_18_16[6],
                   '3600':  adev_18_16[7]}
        #df_results_d18 = df_results_d18.append(new_row, ignore_index = True)
        d18O_rows_list.append(new_row)
 
        # Save data into dataframes dD
        new_row = {'H2O': Picarro_data_subset['H2O'].mean(),
                   'std':Picarro_data_subset['H2O'].std(),
                   'run lenght':(Picarro_data_subset.count()[0]/curr_rate)/3600,
                   '1':     adev_D_H[0],
                   '10':    adev_D_H[1],
                   '60':    adev_D_H[2],
                   '120':   adev_D_H[3],
                   '270':   adev_D_H[4],
                   '600':   adev_D_H[5],
                   '1800':  adev_D_H[6],
                   '3600':  adev_D_H[7]}
        #df_results_dD = df_results_dD.append(new_row, ignore_index = True)
        dD_rows_list.append(new_row)
        
        # Save data into dataframes dx
        new_row = {'H2O': Picarro_data_subset['H2O'].mean(),
                   'std':Picarro_data_subset['H2O'].std(),
                   'run lenght':(Picarro_data_subset.count()[0]/curr_rate)/3600,
                   '1':     adev_dX[0],
                   '10':    adev_dX[1],
                   '60':    adev_dX[2],
                   '120':   adev_dX[3],
                   '270':   adev_dX[4],
                   '600':   adev_dX[5],
                   '1800':  adev_dX[6],
                   '3600':  adev_dX[7]}
        #df_results_dx = df_results_dx.append(new_row, ignore_index = True)
        dx_rows_list.append(new_row)
        
        # Save data into dataframes 17X
        new_row = {'H2O': Picarro_data_subset['H2O'].mean(),
                   'std':Picarro_data_subset['H2O'].std(),
                   'run lenght':(Picarro_data_subset.count()[0]/curr_rate)/3600,
                   '1':     adev_17X_calc[0],
                   '10':    adev_17X_calc[1],
                   '60':    adev_17X_calc[2],
                   '120':   adev_17X_calc[3],
                   '270':   adev_17X_calc[4],
                   '600':   adev_17X_calc[5],
                   '1800':  adev_17X_calc[6],
                   '3600':  adev_17X_calc[7]}
        #df_results_17X = df_results_17X.append(new_row, ignore_index = True)
        D17X_calc_rows_list.append(new_row)
        
        # Save ramp data
        new_row = {'ID':np.repeat(metadata['ID'][sample], Picarro_data_subset['H2O'].count()),
                   'H2O' : Picarro_data_subset['H2O'].to_numpy()}
        
        df_ramp_H2O = pd.concat([df_ramp_H2O, pd.DataFrame(new_row)])
        
        # Save linear regression paramters
        reg_mod = linregress(np.arange(0, len(Picarro_data_subset.index)), Picarro_data_subset['H2O'])
        buf = [Picarro_data_subset['H2O'].mean(), reg_mod.slope, reg_mod.pvalue]
        lin_regress_params = np.append(lin_regress_params, buf)

# Drop first row of dataframes
#df_results_H2O.drop(index = 0, inplace = True)
#df_results_d17.drop(index = 0, inplace = True)
#df_results_d18.drop(index = 0, inplace = True)
#df_results_dD.drop(index = 0, inplace = True)
#df_results_17X.drop(index = 0, inplace = True)
#df_results_dx.drop(index = 0, inplace = True)
# Create dataframes
df_results_H2O  = pd.DataFrame(H2O_rows_list)
df_results_d17  = pd.DataFrame(d17O_rows_list)
df_results_d18  = pd.DataFrame(d18O_rows_list)
df_results_dD   = pd.DataFrame(dD_rows_list)
df_results_dx   = pd.DataFrame(dx_rows_list)
df_results_17X  = pd.DataFrame(D17X_calc_rows_list)

# Reshape linear regression paramters table
lin_regress_params = np.reshape(lin_regress_params, (int(np.floor(len(lin_regress_params)/3)), 3))
lin_regress_params[:,1] = lin_regress_params[:,1]*3600

# Create dummy index for H2O ramp
ts = np.diff(Picarro_data_subset.index).mean()/np.timedelta64(1, 's')
dummy_index = np.arange(0, ts*df_ramp_H2O.count()[0], ts)
df_ramp_H2O.index = dummy_index
#drop first entry of df
df_ramp_H2O.drop(index = 0, axis = 1, inplace = True)

#%% Plot 17X as a function of humidity
df_results_17X.sort_values('H2O', inplace = True)
RSD = df_results_H2O['std']/df_results_H2O['H2O']
mask = RSD<0.2
fig, ax = plt.subplots(dpi = 300)
ax.plot(df_results_17X['H2O'][mask], df_results_17X['600'][mask], 'k-o', label = 't = 600 s')
ax.plot(df_results_17X['H2O'][mask], df_results_17X['1800'][mask], 'k--o', label = 't = 1800 s')
ax.plot(df_results_17X['H2O'][mask], df_results_17X['3600'][mask], 'k:o', label = 't = 3600 s')
ax.set_xlabel('H$_2$O (ppmv)')
ax.set_ylabel('$\sigma_{Allan}$ $\Delta^{17}$O (per meg)')
ax.legend()
ax.grid(which = 'minor')
ax.set_ylim([0,20])
plt.tight_layout()

#%% Plot ramp
props = dict(boxstyle='round', facecolor='lavender', alpha=0.3)
fig, ax = plt.subplots(dpi = 300, figsize = (4, 4*0.8069336521219367))
divider = 3600
ax.plot(df_ramp_H2O.index/divider, df_ramp_H2O['H2O'], 'k')
ax.grid(which = 'minor')
for ID in metadata['ID']:
    avg = round(np.mean(df_ramp_H2O['H2O'][df_ramp_H2O['ID'] == ID]))
    std = round(np.std(df_ramp_H2O['H2O'][df_ramp_H2O['ID'] == ID]))
    #buff_str = "%d ± %d (%.1f%%)" % (avg, std, 100*std/avg)
    buff_str = "%d ± %d" % (avg, std)
    # Interval
    xvals = [df_ramp_H2O.index[df_ramp_H2O['ID'] == ID][0], 
             df_ramp_H2O.index[df_ramp_H2O['ID'] == ID][-1]]
    ax.plot(np.array(xvals)/divider, [avg, avg], 'r-', lw = 1.5, alpha = 0.75)
    xval = xvals[1]/divider+1
    if ID == "BER11k":
        yval = avg - 300
    elif ID == "BER11kbis":
        yval = avg + 300
    elif ID == "BER0k5":
        xval = (xvals[1]-1000)/divider+1
        yval = avg - 400
    else:
        yval = avg
        
    ax.text(xval, yval, buff_str, color = 'r', size = 8)#, bbox = props)

ax.set_xticks(np.arange(-20, 100, 20))
ax.set_xticklabels(['','0', '20', '40', '60', '80'])
ax.set_xlabel('Time (hours)')
ax.set_ylabel('H$_{2}$O (ppmv)')

#axins = inset_axes(ax, width="33%", height="33%", loc=3, borderpad=2)
axins = ax.inset_axes([0.11, 0.09, 0.40, 0.40])


axis_y_lim = [-35, 35]

axins.plot(lin_regress_params[:,0], lin_regress_params[:,1],
           'k-o', lw = 1,
           alpha = 1)
axins.set_ylim(axis_y_lim)
axins.set_xlim([-200, 21000])
#axins.set_xticks([500, 10000, 20000], fontsize=3, rotation=90)
axins.grid(which = 'minor')
axins.set_facecolor('whitesmoke')
axins.text(550, -30, "$\Delta H_{2}O/\Delta T$ slope\n (ppmv/h)", fontsize = 8)
plt.tight_layout()
fig.savefig('asd.pdf')
