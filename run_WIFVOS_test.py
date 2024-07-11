#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:56:09 2022

A versatile water vapor generation module for vapor isotope
calibration and liquid isotope measurements

by Steen-Larsen and Zannoni

Atmospheric Measurement Techniques

Code for reproducing Figure 7

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

from import_picarro_raw_data import import_picarro_raw_data

from lmfit import Model
from lmfit.lineshapes import exponential, powerlaw
from lmfit.models import ExpressionModel, PolynomialModel


plot_in_local_time = True
local_offset = +1 # UTC + 2

#%%Style
from matplotlib import style as mstyle
mstyle.use('article.mplstyle')

colors = {'H2O'         :'black',
          'd17O'        :'dodgerblue',
          'd18O'        :'forestgreen',
          'dD'          :'chocolate',
          'd-excess'    :'black',
          'O17Excess'   :'black',}

BER_col = 'r'
FIN_col = 'b'
GLW_col = 'k'
alpha_bars = .99

#%% Load Picarro data
paths = ['../DATA/HKDS2092/2022/04/07/',
         '../DATA/HKDS2092/2022/04/10/',
         '../DATA/HKDS2092/2022/04/13/',
         '../DATA/HKDS2092/2022/04/14/',
         ]# Directory for daily data


#%% List files and import
for file_path in paths:
    if 'Picarro_data' in locals():
        new_Picarro_data = import_picarro_raw_data(file_path, '17O')
        Picarro_data = pd.concat([Picarro_data, new_Picarro_data])
        del(new_Picarro_data)
    else:
        Picarro_data = import_picarro_raw_data(file_path, '17O')

#%% Calibrate humidity
# Factors calcualted with VAISALA HMP113 
slope       = 0.733863
intercept   = 33.805885
Picarro_data['H2O'] = Picarro_data['H2O']*slope + intercept

#%% Load metadata for analysis
metadata_file_name = 'Timings_Picarro.xlsx'
metadata = pd.read_excel(io=metadata_file_name, sheet_name='WIFVOS_test')
print(metadata.head(5))  # print first 5 rows of the dataframe

#%% Create dataframe to store values
results = pd.DataFrame(columns = ['ID', 'Date', 'H2O', 'H2Oerr', 'd17O', 'd17Oerr', 'd18O', 'd18Oerr', 'dD', 'dDerr', 'n', 'STD'])

#%% Run analysis per line
for idx, sample in zip(range(len(metadata.index)), metadata.index):
    start_winOI = metadata['Avg_start'][sample]
    stop_winOI  = metadata['Avg_stop'][sample]
    Picarro_date_mask = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)
    Picarro_data_subset = Picarro_data[Picarro_date_mask]
    # Append data to dataframe        
    new_row = {'ID':        metadata['ID'][sample],
               'Date':      metadata['Avg_start'][sample], 
               'H2O':       Picarro_data_subset['H2O'].mean(), 
               'H2Oerr':    Picarro_data_subset['H2O'].std(),
               'H2Osem':    Picarro_data_subset['H2O'].sem(),
               'd17O':      Picarro_data_subset['Delta_17_16'].mean(),
               'd17Oerr':   Picarro_data_subset['Delta_17_16'].std(),
               'd17Osem':   Picarro_data_subset['Delta_17_16'].sem(),
               'd18O':      Picarro_data_subset['Delta_18_16'].mean(),
               'd18Oerr':   Picarro_data_subset['Delta_18_16'].std(),
               'd18Osem':   Picarro_data_subset['Delta_18_16'].sem(),
               'dD':        Picarro_data_subset['Delta_D_H'].mean(),
               'dDerr':     Picarro_data_subset['Delta_D_H'].std(),
               'dDsem':     Picarro_data_subset['Delta_D_H'].sem(),
               'n':         Picarro_data_subset['Delta_17_16'].count(),
               'STD':       metadata['Notes'][sample]}
    #results = results.append(new_row, ignore_index = True)
    results = pd.concat([results, pd.DataFrame(data=new_row, index=[idx])], ignore_index=True)

#%%
mod = ExpressionModel("off + a*exp(b*x)", independent_vars=['x'])
x_pred = np.linspace(500, 3500, 100)

#%% Plot SE instead of STD
fig, ((ax_d18O, ax_dD), (ax_d18Ostd, ax_dDstd)) = plt.subplots(nrows = 2, ncols = 2, sharex = True, figsize = (4, 3), 
                                                               dpi = 300)#, gridspec_kw={'height_ratios': [1,2,2]})


mask = (results['STD'] == 'BER') & (results['H2O'] < 3500)
ax_d18O.errorbar(x = results['H2O'][mask], y=results['d18O'][mask].iloc[0] - results['d18O'][mask], 
                 #xerr = results['H2Oerr'][mask]/np.sqrt(len(results['d18Oerr'][mask])),
                 #yerr=results['d18Oerr'][mask]/np.sqrt(len(results['d18Oerr'][mask])), 
                 xerr = results['H2Osem'][mask],
                 yerr = results['d18Osem'][mask],
                 markerfacecolor='w',color=BER_col,capsize=2, linestyle='None', alpha = alpha_bars)
y_val = (results['d18O'][mask].iloc[0] - results['d18O'][mask]).values
x_val = results['H2O'][mask].values
result_d18O_BER = mod.fit(y_val, x = x_val, off = 0, a = -1e-4, b = -1e-4)
y_pred = result_d18O_BER.eval(x = x_pred)
ax_d18O.plot(x_pred, y_pred, color=BER_col)

ax_d18Ostd.errorbar(x = results['H2O'][mask], y=results['d18Oerr'][mask], 
                    xerr = results['H2Oerr'][mask],
                    markerfacecolor='w',color=BER_col,capsize=2, linestyle='None', alpha = alpha_bars)

x_std_data = results['H2O'][mask].values
y_std_data = results['d18Oerr'][mask].values

mask = (results['STD'] == 'FIN') & (results['H2O'] < 3500)
ax_d18O.errorbar(x = results['H2O'][mask], y=results['d18O'][mask].iloc[0] - results['d18O'][mask], 
                 xerr = results['H2Osem'][mask],
                 yerr = results['d18Osem'][mask],
                 markerfacecolor='w',color=FIN_col,capsize=2, linestyle='None', alpha = alpha_bars)
y_val = (results['d18O'][mask].iloc[0] - results['d18O'][mask]).values
x_val = results['H2O'][mask].values
result_d18O_FIN = mod.fit(y_val, x = x_val, off = 0, a = -1e-4, b = -1e-4)
y_pred = result_d18O_FIN.eval(x = x_pred)
ax_d18O.plot(x_pred, y_pred, color=FIN_col)

ax_d18Ostd.errorbar(x = results['H2O'][mask], y=results['d18Oerr'][mask], 
                    xerr = results['H2Oerr'][mask],
                    markerfacecolor='w',color=FIN_col,capsize=2, linestyle='None', alpha = alpha_bars)

x_std_data = np.append(x_std_data, results['H2O'][mask].values)
y_std_data = np.append(y_std_data, results['d18Oerr'][mask].values)

mask = (results['STD'] == 'GLW') & (results['H2O'] < 3500)
ax_d18O.errorbar(x = results['H2O'][mask], y=results['d18O'][mask].iloc[0] - results['d18O'][mask], 
                 xerr = results['H2Osem'][mask],
                 yerr = results['d18Osem'][mask],  
                 markerfacecolor='w',color=GLW_col,capsize=2, linestyle='None', alpha = alpha_bars)
y_val = (results['d18O'][mask].iloc[0] - results['d18O'][mask]).values
x_val = results['H2O'][mask].values
result_d18O_GLW = mod.fit(y_val, x = x_val, off = 0, a = -1e-4, b = -1e-4)
y_pred = result_d18O_GLW.eval(x = x_pred)
ax_d18O.plot(x_pred, y_pred, color=GLW_col)

ax_d18Ostd.errorbar(x = results['H2O'][mask], y=results['d18Oerr'][mask], 
                    xerr = results['H2Oerr'][mask],
                    markerfacecolor='w',color=GLW_col,capsize=2, linestyle='None', alpha = alpha_bars)
x_std_data = np.append(x_std_data, results['H2O'][mask].values)
y_std_data = np.append(y_std_data, results['d18Oerr'][mask].values)
                       
# Generate model
result_d18O_std = mod.fit(y_std_data, x = x_std_data, off = 0.2, a = 1e0, b = -1e-5)
y_pred = result_d18O_std.eval(x = x_pred)
# Plot std model
ax_d18Ostd.plot(x_pred, y_pred, 'k--', lw = 1)
print('d18O STD')
print(result_d18O_std.params['off'])
print(result_d18O_std.params['a'])
print(result_d18O_std.params['b'])
print(1 - result_d18O_std.residual.var() / np.var(y_val))


ax_d18O.set_ylim([-2,2])
ax_d18O.grid(which = 'minor')
#ax_d18O.set_ylabel('$\delta^{18}$O$_{diff}$ (‰)')
ax_d18O.set_ylabel('$\Delta\delta^{18}$O (‰)')

ax_d18Ostd.set_ylim([0, 2.1])
ax_d18Ostd.grid(which = 'minor')
ax_d18Ostd.set_ylabel('$\delta^{18}$O 1$\sigma$ (‰)')
ax_d18Ostd.set_xlabel('H$_2$O (ppmv)')




mask = (results['STD'] == 'BER') & (results['H2O'] < 3500)
ax_dD.errorbar(x = results['H2O'][mask], y=results['dD'][mask].iloc[0] - results['dD'][mask], 
               xerr = results['H2Osem'][mask],
               yerr = results['dDsem'][mask],
               markerfacecolor='w',color=BER_col,capsize=2, linestyle='None', alpha = alpha_bars)
y_val = (results['dD'][mask].iloc[0] - results['dD'][mask]).values
x_val = results['H2O'][mask].values
result_dD_BER = mod.fit(y_val, x = x_val, off = 0, a = 1e-4, b = -1e-4)
y_pred = result_dD_BER.eval(x = x_pred)
ax_dD.plot(x_pred, y_pred, color=BER_col)

ax_dDstd.errorbar(x = results['H2O'][mask], y=results['dDerr'][mask], 
                  xerr = results['H2Oerr'][mask],
                  markerfacecolor='w',color=BER_col,capsize=2, linestyle='None', alpha = alpha_bars)
x_std_data = results['H2O'][mask].values
y_std_data = results['dDerr'][mask].values

mask = (results['STD'] == 'FIN') & (results['H2O'] < 3500)
ax_dD.errorbar(x = results['H2O'][mask], y=results['dD'][mask].iloc[0] - results['dD'][mask], 
               xerr = results['H2Osem'][mask],
               yerr = results['dDsem'][mask],
               markerfacecolor='w',color=FIN_col,capsize=2, linestyle='None', alpha = alpha_bars)
y_val = (results['dD'][mask].iloc[0] - results['dD'][mask]).values
x_val = results['H2O'][mask].values
result_dD_FIN = mod.fit(y_val, x = x_val, off = 0, a = 1e-4, b = -1e-4)
y_pred = result_dD_FIN.eval(x = x_pred)
ax_dD.plot(x_pred, y_pred, color=FIN_col)

ax_dDstd.errorbar(x = results['H2O'][mask], y=results['dDerr'][mask], 
                  xerr = results['H2Oerr'][mask],
                  markerfacecolor='w',color=FIN_col,capsize=2, linestyle='None', alpha = alpha_bars)
x_std_data = np.append(x_std_data, results['H2O'][mask].values)
y_std_data = np.append(y_std_data, results['dDerr'][mask].values)


mask = (results['STD'] == 'GLW') & (results['H2O'] < 3500)
ax_dD.errorbar(x = results['H2O'][mask], y=results['dD'][mask].iloc[0] - results['dD'][mask], 
               xerr = results['H2Osem'][mask],
               yerr = results['dDsem'][mask],
               markerfacecolor='w',color=GLW_col, capsize=2, linestyle='None', alpha = alpha_bars)
y_val = (results['dD'][mask].iloc[0] - results['dD'][mask]).values
x_val = results['H2O'][mask].values
result_dD_GLW = mod.fit(y_val, x = x_val, off = 0, a = 1e-4, b = -1e-4)
y_pred = result_dD_GLW.eval(x = x_pred)
ax_dD.plot(x_pred, y_pred, color=GLW_col)

ax_dDstd.errorbar(x = results['H2O'][mask], y=results['dDerr'][mask], 
                  xerr = results['H2Oerr'][mask],
                  markerfacecolor='w',color=GLW_col,capsize=2, linestyle='None', alpha = alpha_bars)
x_std_data = np.append(x_std_data, results['H2O'][mask].values)
y_std_data = np.append(y_std_data, results['dDerr'][mask].values)

# Generate model
result_dD_std = mod.fit(y_std_data, x = x_std_data, off = 0.2, a = 1e0, b = -1e-5)
y_pred = result_dD_std.eval(x = x_pred)
# Plot std model
ax_dDstd.plot(x_pred, y_pred, 'k--', lw = 1)
print('dD STD')
print(result_dD_std.params['off'])
print(result_dD_std.params['a'])
print(result_dD_std.params['b'])
print(1 - result_dD_std.residual.var() / np.var(y_val))


ax_dD.set_ylim([-7, 7])
ax_dD.grid(which = 'minor')
ax_dD.set_ylabel('$\Delta\delta$D (‰)')
#ax_dD.set_ylabel('$delta$D$_{diff}$ (‰)')

ax_dDstd.set_ylim([0, 8])
ax_dDstd.grid(which = 'minor')
ax_dDstd.set_ylabel('$\delta$D 1$\sigma$ (‰)')
ax_dDstd.set_xlabel('H$_2$O (ppmv)')

plt.tight_layout(h_pad = 0.2)#(pad=0.4, w_pad=0.5, h_pad=1.0)