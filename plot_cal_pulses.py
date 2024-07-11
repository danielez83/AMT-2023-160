#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:56:09 2022

A versatile water vapor generation module for vapor isotope
calibration and liquid isotope measurements

by Steen-Larsen and Zannoni

Atmospheric Measurement Techniques

Code for reproducing Figure 9

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
from scipy.signal import find_peaks

from import_picarro_raw_data import import_picarro_raw_data

#%%Style
from matplotlib import style as mstyle
mstyle.use('article.mplstyle')

from defines import colors
marker_sel      = 'o'
marker_multi    = '^'

#%% Load train pulses
show_selector       = False
show_multioven      = True

df_selector         = pd.read_csv('Cal_Pulses_Selector.csv')
# Drop the first two rows
df_selector.drop(index = [0,1], inplace = True)
df_selector.index   = df_selector['Time']
df_selector['d-excess mean'] = df_selector['dD mean'] - 8*df_selector['d18O mean']

df_multioven    = pd.read_csv('Cal_Pulses_MultiOven_new20230609.csv')
df_multioven.index = df_multioven['Time']

#%% Plot data

xlimits = [0,48]
fig, ax = plt.subplots(nrows=5, ncols = 2, sharex= True, dpi = 300, figsize = (4,4))

# H2O STD1 --------------------------------------------------------------------
# ax[0,0].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
#                  df_selector[df_selector['Standard'] == 'SP']['H2O mean'], 
#                  yerr = df_selector[df_selector['Standard'] == 'SP']['H2O std'], 
#                  marker = marker_sel, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,0].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
#                  df_multioven[df_multioven['Standard'] == 'SP']['H2O mean'], 
#                  yerr = df_multioven[df_multioven['Standard'] == 'SP']['H2O std'], 
#                  marker = marker_multi, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,0].grid(which = 'both')
# ax[0,0].set_ylim([9000, 11000])

# # H2O STD2 --------------------------------------------------------------------
# ax[0,1].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
#                  df_selector[df_selector['Standard'] == 'BER']['H2O mean'], 
#                  yerr = df_selector[df_selector['Standard'] == 'BER']['H2O std'], 
#                  marker = marker_sel, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,1].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
#                  df_multioven[df_multioven['Standard'] == 'BER']['H2O mean'], 
#                  yerr = df_multioven[df_multioven['Standard'] == 'BER']['H2O std'], 
#                  marker = marker_multi, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,1].grid(which = 'both')
#ax[0,1].set_ylim([9000, 11000])

# d17O STD1 --------------------------------------------------------------------
if show_selector:
    ax[0,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['d17O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['d17O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[0,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['d17O mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['d17O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
ax[0,1].grid(which = 'both')
#ax[0,0].set_ylim([-30.45, -30.25])
ax[0,0].set_ylabel('$\delta^{17}$O (‰)', size = 7)
# Panel labels
ax[0,0].text(0.90, 0.15, '(a)', horizontalalignment='center',
             verticalalignment='center', transform=ax[0,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[0,1].text(0.90, 0.15, '(b)', horizontalalignment='center',
             verticalalignment='center', transform=ax[0,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))

# d17O STD2 --------------------------------------------------------------------
if show_selector:
    ax[0,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['d17O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['d17O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[0,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['d17O mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['d17O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
ax[0,0].grid(which = 'both')
#ax[0,1].set_ylim([-0.3,-0.1])


# d18O STD1 --------------------------------------------------------------------
if show_selector:
    ax[1,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['d18O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['d18O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[1,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['d18O mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['d18O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
ax[1,1].grid(which = 'both')
#ax[1,0].set_ylim([-56.9,-56.6])
ax[1,0].set_ylabel('$\delta^{18}$O (‰)', size = 7)
# Panel labels
ax[1,0].text(0.90, 0.15, '(c)', horizontalalignment='center',
             verticalalignment='center', transform=ax[1,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[1,1].text(0.90, 0.15, '(d)', horizontalalignment='center',
             verticalalignment='center', transform=ax[1,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))


# d18O STD2 --------------------------------------------------------------------
if show_selector:
    ax[1,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['d18O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['d18O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[1,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['d18O mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['d18O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
ax[1,0].grid(which = 'both')
#ax[1,1].set_ylim([-1.1,-0.8])

# dD STD1 --------------------------------------------------------------------
if show_selector:
    ax[2,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['dD mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['dD std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[2,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['dD mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['dD std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
ax[2,1].grid(which = 'both')
#ax[2,0].set_ylim([-456,-449])
ax[2,0].set_ylabel('$\delta$D (‰)', size = 7)
# Panel labels
ax[2,0].text(0.90, 0.15, '(e)', horizontalalignment='center',
             verticalalignment='center', transform=ax[2,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[2,1].text(0.90, 0.15, '(f)', horizontalalignment='center',
             verticalalignment='center', transform=ax[2,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))



# dD STD2 --------------------------------------------------------------------
if show_selector:
    ax[2,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['dD mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['dD std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[2,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['dD mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['dD std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
ax[2,0].grid(which = 'both')
#ax[2,1].set_ylim([-5, -4])

# d-Oexcess STD1 --------------------------------------------------------------------
if show_selector:
    ax[3,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['d-excess mean'], 
                     yerr = 0,#df_selector[df_selector['Standard'] == 'SP']['d-excess std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[3,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['d-excess mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['d-excess std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
ax[3,1].grid(which = 'both')
#ax[3,0].set_ylim([-0.35,0.25])
ax[3,0].set_ylabel('d-excess (‰)', size = 7)
# Panel labels
ax[3,0].text(0.90, 0.15, '(g)', horizontalalignment='center',
             verticalalignment='center', transform=ax[3,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[3,1].text(0.90, 0.15, '(h)', horizontalalignment='center',
             verticalalignment='center', transform=ax[3,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))



# d-excess STD2 --------------------------------------------------------------------
if show_selector:
    ax[3,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['d-excess mean'], 
                     yerr = 0,#df_selector[df_selector['Standard'] == 'BER']['d-excess std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[3,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['d-excess mean'], 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['d-excess std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
ax[3,0].grid(which = 'both')
#ax[3,1].set_ylim([-0.1, 0.55])


# 17Oexcess STD1 --------------------------------------------------------------------
if show_selector:
    ax[4,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['Ex17 mean']*1000, 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['Ex17 std']*1000, 
                     marker = marker_sel, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[4,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['Ex17 mean']*1000, 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['Ex17 std']*1000, 
                     marker = marker_multi, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
ax[4,1].grid(which = 'both')
#ax[4,0].set_ylim([-0.35,0.25])
ax[4,0].set_ylabel('$\Delta^{17}$O (per meg)', size = 7)
# Panel labels
ax[4,0].text(0.90, 0.15, '(i)', horizontalalignment='center',
             verticalalignment='center', transform=ax[4,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[4,1].text(0.90, 0.15, '(j)', horizontalalignment='center',
             verticalalignment='center', transform=ax[4,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))


# 17Oexcess STD2 --------------------------------------------------------------------
if show_selector:
    ax[4,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['Ex17 mean']*1000, 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['Ex17 std']*1000, 
                     marker = marker_sel, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[4,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['Ex17 mean']*1000, 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['Ex17 std']*1000, 
                     marker = marker_multi, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
ax[4,0].grid(which = 'both')
#ax[4,1].set_ylim([-0.1, 0.55])



ax[4,1].set_xticks(np.arange(xlimits[0], xlimits[1]+1, 12))
ax[4,1].set_xlim([xlimits[0]-2, xlimits[1]+2])
ax[4,0].set_xlim([xlimits[0]-2, xlimits[1]+2])
ax[4,1].set_xlabel('Time (h)')
ax[4,0].set_xlabel('Time (h)')

ax[0,0].set_title('BER')
ax[0,1].set_title('SP')

plt.tight_layout(h_pad = .1, w_pad = .1)
#fig.savefig('asd.pdf')

#%% Plot data - R1

xlimits = [0,48]
fig, ax = plt.subplots(nrows=5, ncols = 2, sharex= True, dpi = 300, figsize = (4,4))

# H2O STD1 --------------------------------------------------------------------
# ax[0,0].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
#                  df_selector[df_selector['Standard'] == 'SP']['H2O mean'], 
#                  yerr = df_selector[df_selector['Standard'] == 'SP']['H2O std'], 
#                  marker = marker_sel, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,0].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
#                  df_multioven[df_multioven['Standard'] == 'SP']['H2O mean'], 
#                  yerr = df_multioven[df_multioven['Standard'] == 'SP']['H2O std'], 
#                  marker = marker_multi, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,0].grid(which = 'both')
# ax[0,0].set_ylim([9000, 11000])

# # H2O STD2 --------------------------------------------------------------------
# ax[0,1].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
#                  df_selector[df_selector['Standard'] == 'BER']['H2O mean'], 
#                  yerr = df_selector[df_selector['Standard'] == 'BER']['H2O std'], 
#                  marker = marker_sel, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,1].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
#                  df_multioven[df_multioven['Standard'] == 'BER']['H2O mean'], 
#                  yerr = df_multioven[df_multioven['Standard'] == 'BER']['H2O std'], 
#                  marker = marker_multi, 
#                  markerfacecolor='w',
#                  color=colors['H2O'], capsize=2, linestyle='None')#, alpha = 0.5)
# ax[0,1].grid(which = 'both')
#ax[0,1].set_ylim([9000, 11000])

# d17O STD1 --------------------------------------------------------------------
ylim_d17O = [-0.075, 0.075]
if show_selector:
    ax[0,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['d17O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['d17O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[0,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['d17O mean']-df_multioven[df_multioven['Standard'] == 'SP']['d17O mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['d17O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
ax[0,1].grid(which = 'minor')
ax[0,1].set_ylim(ylim_d17O)
ax[0,1].set_yticklabels([])

# d17O STD2 --------------------------------------------------------------------
if show_selector:
    ax[0,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['d17O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['d17O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[0,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['d17O mean']-df_multioven[df_multioven['Standard'] == 'BER']['d17O mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['d17O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d17O'],
                     color='k', capsize=2, linestyle='None')
ax[0,0].grid(which = 'minor')
ax[0,0].set_ylim(ylim_d17O)
ax[0,0].set_ylabel('$\delta^{17}$O (‰)', size = 7)

# Panel labels
ax[0,0].text(0.90, 0.15, '(a)', horizontalalignment='center',
             verticalalignment='center', transform=ax[0,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[0,1].text(0.90, 0.15, '(b)', horizontalalignment='center',
             verticalalignment='center', transform=ax[0,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))


# d18O STD1 --------------------------------------------------------------------
ylim_d18O = [-0.075, 0.075]
if show_selector:
    ax[1,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['d18O mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['d18O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[1,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['d18O mean']-df_multioven[df_multioven['Standard'] == 'SP']['d18O mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['d18O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
ax[1,1].grid(which = 'minor')
ax[1,1].set_ylim(ylim_d18O)
ax[1,1].set_yticklabels([])
ax[1,0].set_ylabel('$\delta^{18}$O (‰)', size = 7)

# d18O STD2 --------------------------------------------------------------------
if show_selector:
    ax[1,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['d18O mean']-df_selector[df_selector['Standard'] == 'BER']['d18O mean'].mean(), 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['d18O std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[1,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['d18O mean']-df_multioven[df_multioven['Standard'] == 'BER']['d18O mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['d18O std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d18O'],
                     color='k', capsize=2, linestyle='None')
ax[1,0].grid(which = 'minor')

# Panel labels
ax[1,0].text(0.90, 0.15, '(c)', horizontalalignment='center',
             verticalalignment='center', transform=ax[1,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[1,1].text(0.90, 0.15, '(d)', horizontalalignment='center',
             verticalalignment='center', transform=ax[1,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))

# dD STD1 --------------------------------------------------------------------
ylim_dD = [-0.3, 0.3]
if show_selector:
    ax[2,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['dD mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['dD std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[2,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['dD mean']-df_multioven[df_multioven['Standard'] == 'SP']['dD mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['dD std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
ax[2,1].grid(which = 'minor')
ax[2,1].set_ylim(ylim_dD)
ax[2,1].set_yticklabels([])
ax[2,0].set_ylabel('$\delta$D (‰)', size = 7)

# dD STD2 --------------------------------------------------------------------
if show_selector:
    ax[2,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['dD mean'], 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['dD std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[2,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['dD mean']-df_multioven[df_multioven['Standard'] == 'BER']['dD mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['dD std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['dD'],
                     color='k', capsize=2, linestyle='None')
ax[2,0].grid(which = 'minor')
ax[2,0].set_ylim(ylim_dD)

# Panel labels
ax[2,0].text(0.90, 0.15, '(e)', horizontalalignment='center',
             verticalalignment='center', transform=ax[2,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[2,1].text(0.90, 0.15, '(f)', horizontalalignment='center',
             verticalalignment='center', transform=ax[2,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))




# d-excess STD1 --------------------------------------------------------------------
ylim_dx = [-0.55, 0.55]
if show_selector:
    ax[3,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['d-excess mean'], 
                     yerr = 0,#df_selector[df_selector['Standard'] == 'SP']['d-excess std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[3,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'SP']['d-excess mean']-df_multioven[df_multioven['Standard'] == 'SP']['d-excess mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['d-excess std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
ax[3,1].grid(which = 'minor')
ax[3,1].set_ylim(ylim_dx)
ax[3,1].set_yticklabels([])
ax[3,0].set_ylabel('d-excess (‰)', size = 7)


# d-excess STD2 --------------------------------------------------------------------
if show_selector:
    ax[3,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['d-excess mean'], 
                     yerr = 0,#df_selector[df_selector['Standard'] == 'BER']['d-excess std'], 
                     marker = marker_sel, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[3,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     df_multioven[df_multioven['Standard'] == 'BER']['d-excess mean']-df_multioven[df_multioven['Standard'] == 'BER']['d-excess mean'].mean(), 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['d-excess std'], 
                     marker = marker_multi, 
                     markerfacecolor=colors['d-excess'],
                     color='k', capsize=2, linestyle='None')
ax[3,0].grid(which = 'minor')
ax[3,0].set_ylim(ylim_dx)

# Panel labels
ax[3,0].text(0.90, 0.15, '(g)', horizontalalignment='center',
             verticalalignment='center', transform=ax[3,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[3,1].text(0.90, 0.15, '(h)', horizontalalignment='center',
             verticalalignment='center', transform=ax[3,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))


# 17Oexcess STD1 --------------------------------------------------------------------
ylim_ox = [-50, 50]
if show_selector:
    ax[4,1].errorbar(df_selector[df_selector['Standard'] == 'SP']['Time'], 
                     df_selector[df_selector['Standard'] == 'SP']['Ex17 mean']*1000, 
                     yerr = df_selector[df_selector['Standard'] == 'SP']['Ex17 std']*1000, 
                     marker = marker_sel, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[4,1].errorbar(df_multioven[df_multioven['Standard'] == 'SP']['Time'], 
                     (df_multioven[df_multioven['Standard'] == 'SP']['Ex17 mean']-df_multioven[df_multioven['Standard'] == 'SP']['Ex17 mean'].mean())*1000, 
                     yerr = df_multioven[df_multioven['Standard'] == 'SP']['Ex17 std']*1000, 
                     marker = marker_multi, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
ax[4,1].grid(which = 'minor')
ax[4,1].set_ylim(ylim_ox)
ax[4,1].set_yticklabels([])
ax[4,0].set_ylabel('$\Delta^{17}$O (per meg)', size = 7)



# 17Oexcess STD2 --------------------------------------------------------------------
if show_selector:
    ax[4,0].errorbar(df_selector[df_selector['Standard'] == 'BER']['Time'], 
                     df_selector[df_selector['Standard'] == 'BER']['Ex17 mean']*1000, 
                     yerr = df_selector[df_selector['Standard'] == 'BER']['Ex17 std']*1000, 
                     marker = marker_sel, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
if show_multioven:
    ax[4,0].errorbar(df_multioven[df_multioven['Standard'] == 'BER']['Time'], 
                     (df_multioven[df_multioven['Standard'] == 'BER']['Ex17 mean']-df_multioven[df_multioven['Standard'] == 'BER']['Ex17 mean'].mean())*1000, 
                     yerr = df_multioven[df_multioven['Standard'] == 'BER']['Ex17 std']*1000, 
                     marker = marker_multi, 
                     markerfacecolor=colors['O17Excess'],
                     color='k', capsize=2, linestyle='None')
ax[4,0].grid(which = 'minor')
ax[4,0].set_ylim(ylim_ox)



ax[4,1].set_xticks(np.arange(xlimits[0], xlimits[1]+1, 12))
ax[4,1].set_xlim([xlimits[0]-2, xlimits[1]+2])
ax[4,0].set_xlim([xlimits[0]-2, xlimits[1]+2])
ax[4,1].set_xlabel('Time (h)')
ax[4,0].set_xlabel('Time (h)')

ax[0,0].set_title('BER')
ax[0,1].set_title('SP')

# Panel labels
ax[4,0].text(0.90, 0.15, '(i)', horizontalalignment='center',
             verticalalignment='center', transform=ax[4,0].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))
ax[4,1].text(0.90, 0.15, '(j)', horizontalalignment='center',
             verticalalignment='center', transform=ax[4,1].transAxes, weight = 1000, fontsize=12,
             bbox=dict(boxstyle='square,pad=0', fc='w', ec='none'))

plt.tight_layout(h_pad = .1, w_pad = .5)
fig.savefig('asd.pdf')

#%% Compare enrichment with drift
#Maximum drift in 24h following the datasheet 0.2‰ (δ18O), 0.2‰ (δ17O), 0.8‰ (δD)
from scipy.stats import linregress
model_BER_d17O = linregress((df_multioven['Time'][df_multioven['Standard'] == 'BER']), df_multioven['d17O mean'][df_multioven['Standard'] == 'BER'])
model_SP_d17O = linregress((df_multioven['Time'][df_multioven['Standard'] == 'SP']), df_multioven['d17O mean'][df_multioven['Standard'] == 'SP'])
print('Drift for d17O:        %.3f o/oo / day,   %.3f o/oo / day (BER, SP)' % (model_BER_d17O.slope*24, model_SP_d17O.slope*24))

model_BER_d18O = linregress((df_multioven['Time'][df_multioven['Standard'] == 'BER']), df_multioven['d18O mean'][df_multioven['Standard'] == 'BER'])
model_SP_d18O = linregress((df_multioven['Time'][df_multioven['Standard'] == 'SP']), df_multioven['d18O mean'][df_multioven['Standard'] == 'SP'])
print('Drift for d18O:        %.3f o/oo / day,   %.3f o/oo / day (BER, SP)' % (model_BER_d18O.slope*24, model_SP_d18O.slope*24))

model_BER_dD = linregress((df_multioven['Time'][df_multioven['Standard'] == 'BER']), df_multioven['dD mean'][df_multioven['Standard'] == 'BER'])
model_SP_dD = linregress((df_multioven['Time'][df_multioven['Standard'] == 'SP']), df_multioven['dD mean'][df_multioven['Standard'] == 'SP'])
print('Drift for dD:          %.3f o/oo / day,   %.3f o/oo / day (BER, SP)' % (model_BER_dD.slope*24, model_SP_dD.slope*24))

model_BER_d = linregress((df_multioven['Time'][df_multioven['Standard'] == 'BER']), df_multioven['d-excess mean'][df_multioven['Standard'] == 'BER'])
model_SP_d = linregress((df_multioven['Time'][df_multioven['Standard'] == 'SP']), df_multioven['d-excess mean'][df_multioven['Standard'] == 'SP'])
print('Drift for d-excess:   %.3f o/oo / day,  %.3f o/oo / day (BER, SP)' % (model_BER_d.slope*24, model_SP_d.slope*24))

model_BER_17O = linregress((df_multioven['Time'][df_multioven['Standard'] == 'BER']), df_multioven['Ex17 mean'][df_multioven['Standard'] == 'BER'])
model_SP_17O = linregress((df_multioven['Time'][df_multioven['Standard'] == 'SP']), df_multioven['Ex17 mean'][df_multioven['Standard'] == 'SP'])
#print('Trend for 17O excess statistically not significant. p-values of slopes are: %.2f and %.2f (BER and SP)' % (model_BER_17O.pvalue, ))
#%% Detrend standards
#df_multioven_corrected = df_multioven.copy()
correction_BER = df_multioven['Time'][df_multioven['Standard'] == 'BER'] * model_BER_d17O.slope + model_BER_d17O.intercept
BER_std = np.std(df_multioven['d17O mean'][df_multioven['Standard'] == 'BER'] - correction_BER)
correction_SP = df_multioven['Time'][df_multioven['Standard'] == 'SP'] * model_SP_d17O.slope + model_SP_d17O.intercept
SP_std = np.std(df_multioven['d17O mean'][df_multioven['Standard'] == 'SP'] - correction_SP)
print('Standard deviations of the mean values for d17O are %.2f o/oo, %.2f o/oo (BER, SP)' % (BER_std, SP_std))

correction_BER = df_multioven['Time'][df_multioven['Standard'] == 'BER'] * model_BER_d18O.slope + model_BER_d18O.intercept
BER_std = np.std(df_multioven['d18O mean'][df_multioven['Standard'] == 'BER'] - correction_BER)
correction_SP = df_multioven['Time'][df_multioven['Standard'] == 'SP'] * model_SP_d18O.slope + model_SP_d18O.intercept
SP_std = np.std(df_multioven['d18O mean'][df_multioven['Standard'] == 'SP'] - correction_SP)
print('Standard deviations of the mean values for d18O are %.2f o/oo, %.2f o/oo (BER, SP)' % (BER_std, SP_std))

correction_BER = df_multioven['Time'][df_multioven['Standard'] == 'BER'] * model_BER_dD.slope + model_BER_dD.intercept
BER_std = np.std(df_multioven['dD mean'][df_multioven['Standard'] == 'BER'] - correction_BER)
correction_SP = df_multioven['Time'][df_multioven['Standard'] == 'SP'] * model_SP_dD.slope + model_SP_dD.intercept
SP_std = np.std(df_multioven['dD mean'][df_multioven['Standard'] == 'SP'] - correction_SP)
print('Standard deviations of the mean values for dD are %.2f o/oo, %.2f o/oo (BER, SP)' % (BER_std, SP_std))

correction_BER = df_multioven['Time'][df_multioven['Standard'] == 'BER'] * model_BER_d.slope + model_BER_d.intercept
BER_std = np.std(df_multioven['d-excess mean'][df_multioven['Standard'] == 'BER'] - correction_BER)
correction_SP = df_multioven['Time'][df_multioven['Standard'] == 'SP'] * model_SP_d.slope + model_SP_d.intercept
SP_std = np.std(df_multioven['d-excess mean'][df_multioven['Standard'] == 'SP'] - correction_SP)
print('Standard deviations of the mean values for d-excess are %.2f o/oo, %.2f o/oo (BER, SP)' % (BER_std, SP_std))

print(np.std(df_multioven['Ex17 mean'][df_multioven['Standard'] == 'BER']))
print(np.std(df_multioven['Ex17 mean'][df_multioven['Standard'] == 'SP']))