"""
Created on Thu Dec 22 13:39 2022
@author: Jiaxi Wu
"""
import numpy as np
import pandas as pd
import os

def color_band(color, mag, slope_range, width):
    num = 0
    for b in np.arange(25, 35, 0.1):
        for m in np.arange(slope_range[0], slope_range[1], 0.1):
            xx1 = (mag - b) / m
            xx2 = (mag - b) / m + width
            judge = (color < xx2) & (color > xx1)
            if np.sum(judge) > num:
                num = np.sum(judge)
                slope, inter = m, b
    return np.array([slope, inter])

info = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'sn_analysis.csv'))
slope, inter = np.zeros(len(info)), np.zeros(len(info))
red_lim, faint_lim = np.zeros(len(info)), np.zeros(len(info))

# calculate band slope and interception, red and green selection lines
slope_range_f606, slope_range_f555 = [-7, -4.95], [-3, -0.95]
red_lim_f555, width = 0.6, 1
for i in range(len(info)):
    data = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'phot_csv/{:s}.csv'.format(info['field'][i])))
    color, mag = data['inst_vega_mag1']-data['inst_vega_mag2'], data['inst_vega_mag2']

    if (info['filter'][i] == 'F555W'):
        band = color_band(color, mag, slope_range_f555, width)
        slope[i], inter[i] = band[0], band[1]
        red_lim[i] = red_lim_f555 + info['A555/606'][i] - info['A814'][i]
    else:
        band = color_band(color, mag, slope_range_f606, width)
        slope[i], inter[i] = band[0], band[1]
        red_lim[i] = 0.3 + info['A555/606'][i] - info['A814'][i]
        
    mini = np.min(mag)//0.1*0.1 + 0.1
    maxi = np.max(mag)//0.1*0.1 + 0.1
    maglist = np.arange(mini, maxi, 0.1)
    err = np.zeros(len(maglist))
    for j in range(len(maglist)):
        judge = (mag>(maglist[j]-0.1)) & (mag<(maglist[j]+0.1))
        subdata = data[judge].reset_index(drop=True)
        if len(subdata) > 0:
            err[j] = np.mean(subdata['mag2_err'])
    mini = np.min(np.abs(err-0.1))
    faint_lim[i] = maglist[np.abs(err-0.1)==mini][0] + 0.5
    if 'NGC4258' in info['field'][i]:
            faint_lim[i] = 27.5

info['red_lim'], info['faint_lim'], info['slope_bc'], info['inter_bc'] = red_lim, faint_lim, slope, inter
info.to_csv(os.path.join(os.path.split(__file__)[0], 'sn_analysis_new.csv'), index=False)