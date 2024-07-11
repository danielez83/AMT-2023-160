    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:56:09 2022

A versatile water vapor generation module for vapor isotope
calibration and liquid isotope measurements

by Steen-Larsen and Zannoni

Atmospheric Measurement Techniques

Code for reproducing Figure 4 and Figure 10

Update the "paths" list with the actual path of PICARRO data

@author: daniele.zannoni@uib.no
"""

#%% Add path to Python
import os
import sys
#%% Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy.stats import norm, kstest, linregress

from import_picarro_raw_data import import_picarro_raw_data
#%%
calibrate_data      = True
d17O_cal_params     = [0.958390535, -0.461053701]
d18O_cal_params     = [0.96478044, -0.375228501]
dD_cal_params       = [0.937748326, -6.009660321]

#%%Style
from matplotlib import style as mstyle
mstyle.use('article.mplstyle')

from defines import colors

#%% Load Picarro data
paths = ['../DATA/HKDS2156/2022/10/10/',
        '../DATA/HKDS2156/2022/10/11/',
        '../DATA/HKDS2156/2022/10/12/',
        '../DATA/HKDS2156/2022/10/13/',
        '../DATA/HKDS2156/2022/10/14/']# Directory for daily data
# Define time region
start_date_str      = '10-10-2022 16:20'
stop_date_str       = '14-10-2022 12:25'

#%% List files and import
for file_path in paths:
    if "Picarro_data" in locals():
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
    

        
#%% Load injection timins
df_injections = pd.read_csv('SP_BER_inj_time.csv')
df_injections['Start'] = pd.to_datetime(df_injections['Start'], dayfirst = True)
df_injections['Stop'] = pd.to_datetime(df_injections['Stop'], dayfirst = True)

#%% Subset date/time of interest
# Subset data 
start_winOI = pd.to_datetime(start_date_str, dayfirst = True)
stop_winOI  = pd.to_datetime(stop_date_str, dayfirst = True)
Picarro_date_mask = (Picarro_data.index > start_winOI) & (Picarro_data.index < stop_winOI)

Picarro_data_subset = Picarro_data[Picarro_date_mask]

Picarro_data_subset['17OExcess_calc'] = 1e6*(np.log(Picarro_data_subset['Delta_17_16']/1000 + 1) - 0.528*np.log(Picarro_data_subset['Delta_18_16']/1000 + 1))

#%% Convert absolute step change to relative step change
low_std_mask = (Picarro_data_subset.index < df_injections['Stop'][0] - np.timedelta64(600, 's')) & (Picarro_data_subset.index > df_injections['Stop'][0] - np.timedelta64(4200, 's'))
high_std_mask = (Picarro_data_subset.index < df_injections['Stop'][4] - np.timedelta64(600, 's')) & (Picarro_data_subset.index > df_injections['Stop'][4] - np.timedelta64(4200, 's'))
# The following one matches the calibration
#high_std_mask = (Picarro_data_subset.index < df_injections['Stop'][1] - np.timedelta64(18000, 's')) & (Picarro_data_subset.index > df_injections['Stop'][1] - np.timedelta64(21600, 's'))
#high_std_mask = (Picarro_data_subset.index < df_injections['Stop'][5] - np.timedelta64(600, 's')) & (Picarro_data_subset.index > df_injections['Stop'][6] - np.timedelta64(4200, 's'))
# Please note that the final 17OX value of BER changes a bit depending on the interval you select. 24 Hours delay to estimate the true value seems to
# to be OK yelding 82 ± 2 permeg. --> df_injections row #5
#high_std_mask = (Picarro_data_subset.index < df_injections['Stop'][13] - np.timedelta64(600, 's')) & (Picarro_data_subset.index > df_injections['Stop'][13] - np.timedelta64(4200, 's'))

dummy_dx    = Picarro_data_subset['Delta_D_H'] - 8* Picarro_data_subset['Delta_18_16']
dummy_dx    = (dummy_dx - dummy_dx[low_std_mask].mean())/(dummy_dx[high_std_mask].mean() - dummy_dx[low_std_mask].mean())
dummy_d17O  = (Picarro_data_subset['Delta_17_16'] - Picarro_data_subset['Delta_17_16'][low_std_mask].mean())/(Picarro_data_subset['Delta_17_16'][high_std_mask].mean() - Picarro_data_subset['Delta_17_16'][low_std_mask].mean())
dummy_d18O  = (Picarro_data_subset['Delta_18_16'] - Picarro_data_subset['Delta_18_16'][low_std_mask].mean())/(Picarro_data_subset['Delta_18_16'][high_std_mask].mean() - Picarro_data_subset['Delta_18_16'][low_std_mask].mean())
dummy_dD    = (Picarro_data_subset['Delta_D_H'] - Picarro_data_subset['Delta_D_H'][low_std_mask].mean())/(Picarro_data_subset['Delta_D_H'][high_std_mask].mean() - Picarro_data_subset['Delta_D_H'][low_std_mask].mean())
#dummy_17X   = (Picarro_data_subset['Excess_17'] - Picarro_data_subset['Excess_17'][low_std_mask].mean())/(Picarro_data_subset['Excess_17'][high_std_mask].mean() - Picarro_data_subset['Excess_17'][low_std_mask].mean())
dummy_17X   = (Picarro_data_subset['17OExcess_calc'] - Picarro_data_subset['17OExcess_calc'][low_std_mask].mean())/(Picarro_data_subset['17OExcess_calc'][high_std_mask].mean() - Picarro_data_subset['17OExcess_calc'][low_std_mask].mean())
low_std_17X_avg = Picarro_data_subset['17OExcess_calc'][low_std_mask].mean()
high_std_17X_avg = Picarro_data_subset['17OExcess_calc'][high_std_mask].mean()
low_std_17X_sem = Picarro_data_subset['17OExcess_calc'][low_std_mask].sem()
high_std_17X_sem = Picarro_data_subset['17OExcess_calc'][high_std_mask].sem()
step_mask   = (Picarro_data_subset.index < df_injections['Stop'][4]) & (Picarro_data_subset.index > df_injections['Start'][0])
basetime    = np.append(np.array([0]), np.cumsum(np.diff(Picarro_data_subset.index)/np.timedelta64(1, 'h'))) 
zero_idx = np.where(np.abs(Picarro_data_subset.index - pd.to_datetime('2022-10-10 20:17:00')) == (np.abs( Picarro_data_subset.index - pd.to_datetime('2022-10-10 20:17:00'))).min())[0][0]
basetime = (basetime - basetime[zero_idx])-0.0125 # 0.0125 is a lag estimated from the plot. The lag is probably due to the moving average OR to the actual time to switch from SP to BER

#%% Read INJ file
df_INJ = pd.read_csv('HKDS2156_IsoWater_20221116_165037.csv', sep=',', skipinitialspace=True)
df_INJ['Time Code'] = pd.to_datetime(df_INJ['Time Code'])
dt_INJ = np.mean(np.diff(df_INJ['Time Code']))/np.timedelta64(1,'h')
dt_INJ_no_dead_time = 0.075 # don't consider dead time
# Remove lines?
to_rem = df_INJ.index[df_INJ['d(17_16)Mean'].isna()][0]
df_INJ.drop([to_rem,1], inplace=True)
N = df_INJ['d(17_16)Mean'].count()
df_INJ.index = np.arange(0,N)
#basetime_INJ = np.append(np.array([0]), np.cumsum(np.diff(df_INJ['Time Code']))/np.timedelta64(1,'h'))

basetime_INJ = np.arange(0, N*dt_INJ, dt_INJ)
basetime_INJ_ndt = np.arange(0, N*dt_INJ_no_dead_time, dt_INJ_no_dead_time)
zero_idx_INJ = 48
basetime_INJ = basetime_INJ - basetime_INJ[zero_idx_INJ]
basetime_INJ_ndt = basetime_INJ_ndt - basetime_INJ_ndt[zero_idx_INJ]

reference_idx_max = len(basetime_INJ)
low_std_mask = np.arange(zero_idx_INJ-6, zero_idx_INJ)
high_std_mask = np.arange(reference_idx_max - 6, reference_idx_max)
dummy_d17O_INJ = (df_INJ['d(17_16)Mean'] - df_INJ['d(17_16)Mean'][low_std_mask].mean())/(df_INJ['d(17_16)Mean'][high_std_mask].mean() - df_INJ['d(17_16)Mean'][low_std_mask].mean())
dummy_d18O_INJ = (df_INJ['d(18_16)Mean'] - df_INJ['d(18_16)Mean'][low_std_mask].mean())/(df_INJ['d(18_16)Mean'][high_std_mask].mean() - df_INJ['d(18_16)Mean'][low_std_mask].mean())
dummy_dD_INJ = (df_INJ['d(D_H)Mean'] - df_INJ['d(D_H)Mean'][low_std_mask].mean())/(df_INJ['d(D_H)Mean'][high_std_mask].mean() - df_INJ['d(D_H)Mean'][low_std_mask].mean())
# #%% Plot time series
# window_AVG = 540

# fig, axes = plt.subplot_mosaic(
#     [["left_col", "top_right"],
#      ["left_col", "mid_right"],
#      ["left_col", "bottom_right"]],
#     figsize = (6,3),
#     gridspec_kw={'width_ratios': [2.5, 1]}
# )

# #fig, ax = plt.subplots(ncols = 2, figsize = (6,3),
# #                       gridspec_kw={'width_ratios': [2.5, 1]})


# axes['left_col'].plot(Picarro_data_subset.index,
#            Picarro_data_subset['H2O'].rolling(window_AVG).mean(),
#            color = colors['H2O'])
# axes['left_col'].set_ylabel('H$_2$O (ppmv)')
# axes['left_col'].set_ylim([15000, 19500])
# axes['left_col'].grid('on', which='both')
# axes['left_col'].set_xticklabels([])
# axes['left_col'].xaxis.set_major_formatter(mdates.DateFormatter('%d-%H'))
# axes['left_col'].set_xlabel('Time (day-hour)')


# top0, bottom0 = axes['left_col'].get_ylim()
# booly = True
# for inj in range(df_injections['Start'].count()):
#     axes['left_col'].text(df_injections['Start'][inj], 18300, df_injections['T'][inj])
#     if booly:
#             if df_injections['T'][inj] == 1:
#                 curr_col = 'dodgerblue'
#             elif df_injections['T'][inj] == 2:
#                 curr_col = 'palegreen'
#             elif df_injections['T'][inj] == 3:
#                 curr_col = 'plum'
#             elif df_injections['T'][inj] == 4:
#                 curr_col = 'salmon'
#             #curr_col = 'dodgerblue'            
#             axes['left_col'].fill_betweenx([bottom0, top0],
#                                            df_injections['Start'][inj], df_injections['Stop'][inj],
#                                            color = curr_col, alpha=0.25)
 
#     #booly = not(booly)

# axis_x_lim = [-0.5, 24]
# axis_y_lim = [0.99, 1.003]
# #axis_y_lim = [-0.5, 1.5]
# axes['top_right'].plot(basetime[step_mask],
#                        dummy_d17O[step_mask].rolling(window_AVG).mean(),
#                        color = colors['d17O'], label = '$\delta^{17}$O')
# axes['top_right'].plot(basetime_INJ,
#                        dummy_d17O_INJ,
#                        color = 'k', linestyle = '--', label = '$\delta^{17}$O Inj.')
# axes['top_right'].set_xlim(axis_x_lim)
# axes['top_right'].set_ylim(axis_y_lim)
# axes['top_right'].grid('on', which='both')

# axes['mid_right'].plot(basetime[step_mask],
#                        dummy_d18O[step_mask].rolling(window_AVG).mean(),
#                        color = colors['d18O'], label = '$\delta^{18}$O')
# axes['mid_right'].plot(basetime_INJ,
#                        dummy_d18O_INJ,
#                        color = 'k', linestyle = '--', label = '$\delta^{18}$O Inj.')
# axes['mid_right'].set_xlim(axis_x_lim)
# axes['mid_right'].set_ylim(axis_y_lim)
# axes['mid_right'].grid('on', which='both')

# axes['bottom_right'].plot(basetime[step_mask],
#                           dummy_dD[step_mask].rolling(window_AVG).mean(),
#                           color = colors['dD'], label = '$\delta$D')
# axes['bottom_right'].plot(basetime_INJ,
#                        dummy_dD_INJ,
#                        color = 'k', linestyle = '--', label = '$\delta$D Inj.')
# axes['bottom_right'].set_xlim(axis_x_lim)
# axes['bottom_right'].set_ylim(axis_y_lim)
# axes['bottom_right'].grid('on', which='both')
# #%% Plot time series
# window_AVG = 270 # Same resolution of the liquid injections

# fig, axes = plt.subplot_mosaic(
#     [["top_panel", "top_panel", "top_panel"],
#      ["bottom_left", "bottom_mid", "bottom_right"]],
#     figsize = (6,6),
#     gridspec_kw={'height_ratios': [1, 2.5]},
#     dpi = 300
# )

# #fig, ax = plt.subplots(ncols = 2, figsize = (6,3),
# #                       gridspec_kw={'width_ratios': [2.5, 1]})


# axes['top_panel'].plot(Picarro_data_subset.index,
#            Picarro_data_subset['H2O'].rolling(window_AVG).mean(),
#            color = colors['H2O'])
# axes['top_panel'].set_ylabel('H$_2$O (ppmv)')
# axes['top_panel'].set_ylim([15000, 20000])
# axes['top_panel'].grid('on', which='both')
# axes['top_panel'].set_xticklabels([])
# axes['top_panel'].xaxis.set_major_formatter(mdates.DateFormatter('%d-%H'))
# axes['top_panel'].set_xlabel('Time (day-hour)')


# top0, bottom0 = axes['top_panel'].get_ylim()
# booly = True
# for inj in range(df_injections['Start'].count()):
#     axes['top_panel'].text(df_injections['Start'][inj], 18300, df_injections['T'][inj])
#     if booly:
#             if df_injections['T'][inj] == 1:
#                 curr_col = [1,1,1]#'dodgerblue'
#             elif df_injections['T'][inj] == 2:
#                 curr_col = [.75, 0.75,.75]#'palegreen'
#             elif df_injections['T'][inj] == 3:
#                 curr_col = [.5,.5,.5]#'plum'
#             elif df_injections['T'][inj] == 4:
#                 curr_col = curr_col = [.25,.25,.25]#'salmon'
#             #curr_col = 'dodgerblue'            
#             axes['top_panel'].fill_betweenx([bottom0, top0],
#                                            df_injections['Start'][inj], df_injections['Stop'][inj],
#                                            color = curr_col, alpha=0.25)
 
#     #booly = not(booly)

# axis_x_lim = [-0.5, 16]
# axis_y_lim = [0.99, 1.002]
# xticks_lower_panel = np.arange(0,18,4)
# # axis_x_lim = [-0.2, 0.2]
# # axis_y_lim = [-0.2, 1]
# # xticks_lower_panel = np.arange(0,2,.2)
# #axis_y_lim = [-0.5, 1.5]
# axes['bottom_left'].plot(basetime[step_mask],
#                        dummy_d17O[step_mask].rolling(window_AVG).mean(),
#                        color = colors['d17O'], label = '$\delta^{17}$O')
# axes['bottom_left'].plot(basetime_INJ,
#                        dummy_d17O_INJ,
#                        color = 'k', linestyle = '--', label = '$\delta^{17}$O Inj.')
# axes['bottom_left'].set_xlim(axis_x_lim)
# axes['bottom_left'].set_ylim(axis_y_lim)
# axes['bottom_left'].grid('on', which='both')
# axes['bottom_left'].set_xlabel('Time (hours)')
# axes['bottom_left'].set_ylabel('Normalized step response')
# axes['bottom_left'].set_xticks(xticks_lower_panel)
# axes['bottom_left'].text(1, 1.001, '$\delta^{17}$O', color = colors['d17O'])

# axes['bottom_mid'].plot(basetime[step_mask],
#                        dummy_d18O[step_mask].rolling(window_AVG).mean(),
#                        color = colors['d18O'], label = '$\delta^{18}$O')
# axes['bottom_mid'].plot(basetime_INJ,
#                        dummy_d18O_INJ,
#                        color = 'k', linestyle = '--', label = '$\delta^{18}$O Inj.')
# axes['bottom_mid'].set_xlim(axis_x_lim)
# axes['bottom_mid'].set_ylim(axis_y_lim)
# axes['bottom_mid'].grid('on', which='both')
# axes['bottom_mid'].set_xlabel('Time (hours)')
# axes['bottom_mid'].set_xticks(xticks_lower_panel)
# axes['bottom_mid'].text(1, 1.001, '$\delta^{18}$O', color = colors['d18O'])

# axes['bottom_right'].plot(basetime[step_mask],
#                           dummy_dD[step_mask].rolling(window_AVG).mean(),
#                           color = colors['dD'], label = '$\delta$D')
# axes['bottom_right'].plot(basetime_INJ,
#                        dummy_dD_INJ,
#                        color = 'k', linestyle = '--', label = '$\delta$D Inj.')
# axes['bottom_right'].set_xlim(axis_x_lim)
# axes['bottom_right'].set_ylim(axis_y_lim)
# axes['bottom_right'].grid('on', which='both')
# axes['bottom_right'].set_xlabel('Time (hours)')
# axes['bottom_right'].set_xticks(xticks_lower_panel)
# axes['bottom_right'].text(1, 1.001, '$\delta$D', color = colors['dD'])

# # Remove low panels labels
# axes['bottom_mid'].set_yticklabels([])
# axes['bottom_right'].set_yticklabels([])
# plt.tight_layout()

#%% Plot time series, only step change
window_AVG = 270 # Same resolution of the liquid injections

fig, axes = plt.subplot_mosaic(
    [["bottom_left", "bottom_mid", "bottom_right"]],
    figsize = (6,3.7),
   # gridspec_kw={'height_ratios': [1, 2.5]},
    dpi = 300
)

#fig, ax = plt.subplots(ncols = 2, figsize = (6,3),
#                       gridspec_kw={'width_ratios': [2.5, 1]})


# axis_x_lim = [-0.5, 30]
# axis_y_lim = [0.99, 1.002]
# xticks_lower_panel = np.arange(0,32,8)
axis_x_lim = [0, 3]
#xticks_lower_panel = np.linspace(axis_x_lim[0],axis_x_lim[1], 6)
xticks_lower_panel = np.arange(axis_x_lim[0], axis_x_lim[1], 0.5)
axis_y_lim = [0.9, 1.01]
axes['bottom_left'].plot(basetime[step_mask],
                       dummy_d17O[step_mask].rolling(window_AVG).mean(),
                       color = colors['d17O'], #label = '$\delta^{17}$O Cal. Sys.',
                       label = 'Calibration Syst.',
                       alpha = 0.80)
axes['bottom_left'].scatter(basetime_INJ,
                            dummy_d17O_INJ,
                            edgecolor = 'k', facecolor = colors['d17O'], 
                            marker = 'o', s = 8,
                            linewidth = .5,
                            #label = '$\delta^{17}$O Inj.')
                            label = 'Injections')
axes['bottom_left'].scatter(basetime_INJ_ndt,
                            dummy_d17O_INJ,
                            color = colors['d17O'], 
                            marker = 'x', s = 8,
                            linewidth = .5,
                            #label = '$\delta^{17}$O Inj. (w/o d.t)')
                            label = 'Injections (w/o d.t)')
axes['bottom_left'].set_xlim(axis_x_lim)
axes['bottom_left'].set_ylim(axis_y_lim)
axes['bottom_left'].grid('on', which='major')
#axes['bottom_left'].set_xlabel('Time (hours)')
axes['bottom_left'].set_ylabel('Normalized step change')
axes['bottom_left'].set_xticks(xticks_lower_panel)
axes['bottom_left'].minorticks_off()
axes['bottom_left'].legend(fontsize = 8, loc = 'lower right')
axes['bottom_left'].set_title('$\delta^{17}$O')
#axes['bottom_left'].text(1, 1.001, '$\delta^{17}$O', color = colors['d17O'])
# Panel label
axes['bottom_left'].text(.1, .95, '(a)', horizontalalignment='center',
                         verticalalignment='center', transform=axes['bottom_left'].transAxes, weight = 1000, fontsize=14)

axes['bottom_mid'].plot(basetime[step_mask],
                       dummy_d18O[step_mask].rolling(window_AVG).mean(),
                       color = colors['d18O'], #label = '$\delta^{18}$O Cal. Sys.',
                       label = 'Calibration Syst.',
                       alpha = 0.8)
axes['bottom_mid'].scatter(basetime_INJ,
                           dummy_d18O_INJ,
                           edgecolor = 'k', facecolor = colors['d18O'], 
                           marker = 'o', s = 8,
                           linewidth = .5,
                           #label = '$\delta^{18}$O Inj.'
                           label = 'Injections')
axes['bottom_mid'].scatter(basetime_INJ_ndt,
                            dummy_d18O_INJ,
                            color = colors['d18O'], 
                            marker = 'x', s = 8,
                            linewidth = .5,
                            #label = '$\delta^{18}$O Inj. (w/o d.t)'
                            label = 'Injections (w/o d.t)')

axes['bottom_mid'].set_xlim(axis_x_lim)
axes['bottom_mid'].set_ylim(axis_y_lim)
axes['bottom_mid'].grid('on', which='major')
axes['bottom_mid'].set_xlabel('Time (hours)')
axes['bottom_mid'].set_xticks(xticks_lower_panel)
axes['bottom_mid'].minorticks_off()
axes['bottom_mid'].legend(fontsize = 8, loc = 'lower right')
axes['bottom_mid'].set_title('$\delta^{18}$O')
#axes['bottom_mid'].text(1, 1.001, '$\delta^{18}$O', color = colors['d18O'])
# Panel label
axes['bottom_mid'].text(.1, .95, '(b)', horizontalalignment='center',
                         verticalalignment='center', transform=axes['bottom_mid'].transAxes, weight = 1000, fontsize=14)

axes['bottom_right'].plot(basetime[step_mask],
                          dummy_dD[step_mask].rolling(window_AVG).mean(),
                          color = colors['dD'], #label = '$\delta$D Cal. Sys.',
                          label = 'Calibration Syst.',
                          alpha = 0.8)
axes['bottom_right'].scatter(basetime_INJ,
                             dummy_dD_INJ,
                             edgecolor = 'k', facecolor = colors['dD'], 
                             marker = 'o', s = 8,
                             linewidth = .5,
                             label = 'Injections')
axes['bottom_right'].scatter(basetime_INJ_ndt,
                            dummy_dD_INJ,
                            color = colors['dD'], 
                            marker = 'x', s = 8,
                            linewidth = .5,
                            label = 'Injections (w/o d.t)')

axes['bottom_right'].set_xlim(axis_x_lim)
axes['bottom_right'].set_ylim(axis_y_lim)
axes['bottom_right'].grid('on', which='major')
#axes['bottom_right'].set_xlabel('Time (hours)')
axes['bottom_right'].set_xticks(xticks_lower_panel)
axes['bottom_right'].minorticks_off()
axes['bottom_right'].legend(fontsize = 8, loc = 'lower right')
axes['bottom_right'].set_title('$\delta$D')
#axes['bottom_right'].text(1, 1.001, '$\delta$D', color = colors['dD'])
# Panel label
axes['bottom_right'].text(.1, .95, '(c)', horizontalalignment='center',
                         verticalalignment='center', transform=axes['bottom_right'].transAxes, weight = 1000, fontsize=14)

# Remove low panels labels
axes['bottom_mid'].set_yticklabels([])
axes['bottom_right'].set_yticklabels([])
#plt.tight_layout(w_pad = 0.00001)
plt.subplots_adjust(left=0.1,
                    bottom=0.15, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.1, 
                    hspace=0.1)
#plt.tight_layout()
#fig.savefig('asd.pdf')
#%% Print 17O Excess step change
# for this specific plot use a different style
# works only with %matplotlin qt5
mstyle.use('article_2.mplstyle')
window_AVG = 900
window_AVG_samples = 3600
timebase = np.mean(np.diff(dummy_17X[step_mask].index))/np.timedelta64(1, 's')
win_lenght = round(window_AVG_samples/timebase)
# calculate averaged samples
n_samples = int(np.floor(len(dummy_17X)/win_lenght))
samples = np.zeros(n_samples)
x_samples = np.zeros(n_samples)
for sample in range(n_samples):
    idx_start = sample*win_lenght
    samples[sample] = np.mean(dummy_17X[step_mask][idx_start:idx_start+win_lenght])
    x_samples[sample] = np.mean(basetime[step_mask][idx_start:idx_start+win_lenght])
    

#fig, ax = plt.subplots()
#fig, ax = plt.subplots(dpi = 300, figsize = (4, 4*0.8069336521219367))
fig, ax = plt.subplots(dpi = 300, figsize = (3.5, 5*0.8069336521219367))

ax.plot(basetime[step_mask],
        dummy_17X[step_mask].rolling(window_AVG).mean(),
        color = colors['O17Excess'], 
        label = '$\Delta&{17}$O Cal. Sys.',
        alpha = 0.8)

str_buff = "(%d ± %d permeg)" % (low_std_17X_avg, low_std_17X_sem)
#ax.text(-5, -0.23, str_buff, fontsize = 10)
str_buff = "(%d ± %d permeg)" % (high_std_17X_avg, high_std_17X_sem)
#ax.text(1, 1.15, str_buff, fontsize = 10)

axins = inset_axes(ax, width="50%", height="50%", loc=4, borderpad=2.5)

axis_y_lim = [0.85, 1.15]
axins.plot(basetime[step_mask],
           dummy_17X[step_mask].rolling(window_AVG).mean(),
           color = colors['O17Excess'], 
           label = '$\Delta&{17}$O Cal. Sys.',
           alpha = 0.5)
axins.scatter(x_samples, samples, fc = [1,1,1,1], ec = 'k')
axins.set_ylim(axis_y_lim)
axins.set_xlim([-1,24])
axins.set_facecolor('whitesmoke')
axins.grid(which = 'minor')
plt.tight_layout()


ax.set_ylabel('Normalized step change')
ax.set_xlabel('Time (hours)')
ax.grid(which = 'minor')

# Remove nans
#mask_new = ~np.isnan(dummy_17X[step_mask]) & ~np.isnan(Picarro_data_subset['17OExcess_calc'][step_mask])
# add secondary axis (09/01/2024)
def norm_to_permeg(x):
    n_2_p = linregress((dummy_17X[step_mask], Picarro_data_subset['17OExcess_calc'][step_mask]))
    return x*n_2_p.slope+n_2_p.intercept
def permeg_to_norm(x):
    p_2_n = linregress((Picarro_data_subset['17OExcess_calc'][step_mask], dummy_17X[step_mask]))
    return x*p_2_n.slope+p_2_n.intercept
    
secax_y = ax.secondary_yaxis(
    'right', functions=(norm_to_permeg, permeg_to_norm))#, yticks )
secax_y.set_ylabel('$\Delta^{17}$O (permeg)')

secax_y_inset = axins.secondary_yaxis(
    'right', functions=(norm_to_permeg, permeg_to_norm))#, yticks )

plot_cal_int = False
if plot_cal_int:
    std1 = (basetime[step_mask] > -1) & (basetime[step_mask] < 0)
    ax.plot(basetime[step_mask][std1],
            dummy_17X[step_mask][std1].rolling(window_AVG).mean(), 'r-')
    std2 = (basetime[step_mask] > 2) & (basetime[step_mask] < 3)
    ax.plot(basetime[step_mask][std2],
            dummy_17X[step_mask][std2].rolling(window_AVG).mean(), 'r-')
plt.tight_layout()


fig.savefig('asd.pdf')

# #%% Print d-excess step change



# #window_AVG = 3600 
# window_AVG = 300
# window_AVG_samples = 1800
# timebase = np.mean(np.diff(dummy_dx[step_mask].index))/np.timedelta64(1, 's')
# win_lenght = round(window_AVG_samples/timebase)
# # calculate averaged samples
# n_samples = int(np.floor(len(dummy_dx)/win_lenght))
# samples = np.zeros(n_samples)
# x_samples = np.zeros(n_samples)
# for sample in range(n_samples):
#     idx_start = sample*win_lenght
#     samples[sample] = np.mean(dummy_dx[step_mask][idx_start:idx_start+win_lenght])
#     x_samples[sample] = np.mean(basetime[step_mask][idx_start:idx_start+win_lenght])
    

# fig, ax = plt.subplots()
# ax.plot(basetime[step_mask],
#         dummy_dx[step_mask].rolling(window_AVG).mean(),
#         color = colors['d-excess'], 
#         label = 'd-excess Cal. Sys.',
#         alpha = 0.8)

# axins = inset_axes(ax, width="50%", height="50%", loc=4, borderpad=2)

# ax.set_xlim([-.75, 30])
# axis_y_lim = [0.90, 1.1]

# axins.plot(basetime[step_mask],
#            dummy_dx[step_mask].rolling(window_AVG).mean(),
#            color = colors['d-excess'], 
#            label = 'd-excess Cal. Sys.',
#            alpha = 0.5)
# axins.scatter(x_samples, samples, fc = [1,1,1,1], ec = 'k')
# axins.set_ylim(axis_y_lim)
# axins.set_xlim([0, 12])


# ax.set_ylabel('Normalized step change')
# ax.set_xlabel('Time (hours)')
# plt.tight_layout()

#%% Print statistics
for inj in range(df_injections['Start'].count()):           
    mask = (Picarro_data_subset.index < df_injections['Stop'][inj]) & (Picarro_data_subset.index > df_injections['Start'][inj])
    mean_val = Picarro_data_subset['H2O'][mask].mean()
    sigma_val = Picarro_data_subset['H2O'][mask].std()
    n = Picarro_data_subset['H2O'][mask].count()
    dur = (df_injections['Stop'][inj] - df_injections['Start'][inj])/np.timedelta64(3600, 's')
    print('%s, %s, %d, %d, %d, %d' % (df_injections['Start'][inj],
                              df_injections['Stop'][inj],
                              dur,
                              mean_val,
                              sigma_val,
                              n))

#%% Print statistics
print("Statistics -----------------------------------------")
print("H2O:     %d ± %d" % (Picarro_data_subset['H2O'].mean(), Picarro_data_subset['H2O'].std()))
print("d18O:    %.2f ± %.2f" % (Picarro_data_subset['Delta_18_16'].mean(), Picarro_data_subset['Delta_18_16'].std()))
print("dD:      %.2f ± %.2f" % (Picarro_data_subset['Delta_D_H'].mean(), Picarro_data_subset['Delta_D_H'].std()))

#%% Print statistics
levels = [0.63, 0.95, 0.99, 0.995, 0.9999]
print("Timings for continuous injection ----------------------------")
for level in levels:
    print("Level: %.4f" % level)
    t = np.where(dummy_d17O.rolling(window_AVG).mean()>level)[0][0]
    print("d17O: %.4f hours" % basetime[t])
    t = np.where(dummy_d18O.rolling(window_AVG).mean()>level)[0][0]
    print("d18O: %.4f hours" % basetime[t])
    t = np.where(dummy_dD.rolling(window_AVG).mean()>level)[0][0]
    print("dDO: %.4f hours" % basetime[t])
    print('----')
    
print("Timings for discrete injection ------------------------------")
for level in levels:
    print("Level: %.4f" % level)
    t = np.where(dummy_d17O_INJ >level)[0][0]
    print("d17O: %.4f hours" % basetime_INJ[t])
    t = np.where(dummy_d18O_INJ >level)[0][0]
    print("d18O: %.4f hours" % basetime_INJ[t])
    t = np.where(dummy_dD_INJ >level)[0][0]
    print("dDO: %.4f hours" % basetime_INJ[t])
    print('----')

print("Timings for discrete injection W/O DEAD TIME ------------------------------")
for level in levels:
    print("Level: %.4f" % level)
    t = np.where(dummy_d17O_INJ >level)[0][0]
    print("d17O: %.4f hours" % basetime_INJ_ndt[t])
    t = np.where(dummy_d18O_INJ >level)[0][0]
    print("d18O: %.4f hours" % basetime_INJ_ndt[t])
    t = np.where(dummy_dD_INJ >level)[0][0]
    print("dDO: %.4f hours" % basetime_INJ_ndt[t])
    print('----')    
    
#%%
if calibrate_data:
    print('Calibrated data...')