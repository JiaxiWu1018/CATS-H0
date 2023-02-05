"""
Created on Sun Jan 29 22:45 2023
@author: Jiaxi Wu
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import linalg
import os

# first part - spatial clipping
def blue_selection(xx, yy, xmin, ymax):
    judge_blue = (xx <= xmin) & (yy <= ymax)
    judge_red = [not i for i in judge_blue]
    return judge_red, judge_blue

def num_den(xx, yy, bins = 50):
    fig, ax = plt.subplots(figsize = (6, 6))
    hist, binx, biny, imag = ax.hist2d(xx, yy, bins = bins, cmap = cm.binary)
    plt.close()
    return hist, binx, biny

def gloess2d(xx, yy, zz, tau = 0.2):
    yest = np.zeros(zz.shape)
    xlen, ylen = np.arange(xx.shape[0])/10, np.arange(xx.shape[1])/10
    ylen, xlen = np.meshgrid(ylen, xlen)
    w = np.array([])
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            w_new = np.array(np.exp(-((xlen - i/10) ** 2 + (ylen - j/10) ** 2) / (2 * tau ** 2)))
            w = np.append(w, w_new)
    w = w.reshape(xx.shape[0], xx.shape[1], xx.shape[0], xx.shape[1])
    for i in range(zz.shape[0]):
        for j in range(zz.shape[1]):
            weight = w[i, j, :, :]
            A = np.array([[np.sum(weight*xx**2),np.sum(weight*xx*yy),np.sum(weight*xx)],
                          [np.sum(weight*xx*yy),np.sum(weight*yy**2),np.sum(weight*yy)],
                          [np.sum(weight*xx),np.sum(weight*yy),np.sum(weight)]])
            b = np.array([np.sum(weight*xx*zz), np.sum(weight*yy*zz), np.sum(weight*zz)])
            theta = linalg.solve(A, b)
            yest[i, j] = theta[0]*xx[i, j]+theta[1]*yy[i, j]+theta[2]
    return yest

def spatial_clip(den, binx, biny, data, spclip):
    if (spclip == 'no'):
        return data
    if (spclip[-1] == 'p'):
        trun = int(spclip[:-1]) / 100 * max(4, np.max(den))
    else:
        trun = float(spclip)
    print(trun)
    judget = np.array([True]*len(data))
    for i in range(den.shape[0]):
        for j in range(den.shape[1]):
            if den[i, j] >= trun:
                judge = (data['x'] >= binx[i]) & (data['x'] <= binx[i + 1]) & (data['y'] >= biny[j]) & (data['y'] <= biny[j + 1])
                judge = [not item for item in judge]
                judget = judget & judge
    data_clipped = data[judget]
    judget = [not item for item in judget]
    data_cut = data[judget]
    return data_clipped, data_cut
'''
# spatial clipping
def clipping(spclip):
    info = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'sn_analysis_new.csv'))
    for i in range(len(info)):
        print(info['field'][i])
        data = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'phot_csv/{:s}.csv'.format(info['field'][i])))
        if (info['field'][i] == 'NGC1404'):
            data = data[data['y'] > 3400].reset_index(drop=True)
        color, mag = data['inst_vega_mag1']-data['inst_vega_mag2'], data['inst_vega_mag2']
        
        judge_red, judge_blue = blue_selection(color, mag, info['red_lim'][i], info['faint_lim'][i])
        data_red, data_blue = data[judge_red].reset_index(drop=True), data[judge_blue].reset_index(drop=True)
        Z_blue, binx, biny = num_den(np.array(data_blue['x']), np.array(data_blue['y']))
        Z_red, binx, biny = num_den(np.array(data_red['x']), np.array(data_red['y']), (binx, biny))
        X, Y = (binx[1:]+binx[:-1])/2, (biny[1:]+biny[:-1])/2
        Y, X = np.meshgrid(Y, X)
        
        Z_ratio = np.zeros(np.shape(Z_red))
        for j in range(len(Z_blue)):
            for k in range(len(Z_blue[0, :])):
                if Z_blue[j, k] == 0:
                    Z_ratio[j, k] = np.max(Z_red)
                else:
                    Z_ratio[j, k] = Z_red[j, k] / Z_blue[j, k]
        
        Z_smooth = gloess2d(X, Y, Z_blue)
        data_clip, data_cut = spatial_clip(Z_smooth, binx, biny, data, spclip)
        data_clip, data_cut = data_clip.reset_index(drop=True), data_cut.reset_index(drop=True)
        data_clip.to_csv(os.path.join(os.path.split(__file__)[0], 'clipped_csv/{:s}_{:s}.csv'.format(info['field'][i], spclip)), index=False)
        data_cut.to_csv(os.path.join(os.path.split(__file__)[0], 'cut_csv/{:s}_{:s}.csv'.format(info['field'][i], spclip)), index=False)

        fig = plt.figure(figsize = (32, 20))
        ax1 = fig.add_subplot(2, 17, (1,5))
        cs = ax1.contourf(X, Y, Z_ratio, cmap = 'Reds')
        ax1.set_title('Red to Blue Ratio')
        cbar = fig.colorbar(cs)
        ax2 = fig.add_subplot(2, 17, (7,11))
        cs = ax2.contourf(X, Y, Z_blue, cmap = 'Blues')
        ax2.set_title('Blue Star Pop')
        cbar = fig.colorbar(cs)
        ax3 = fig.add_subplot(2, 17, (13,17))
        cs = ax3.contourf(X, Y, Z_smooth, cmap = 'Blues')
        ax3.set_title('Blue Star Pop Smoothed')
        cbar = fig.colorbar(cs)
        
        ax4 = fig.add_subplot(2, 18, (19, 22))
        ax4.plot(data_red['inst_vega_mag1']-data_red['inst_vega_mag2'], data_red['inst_vega_mag2'], 'k.', markersize = 2)
        ax4.plot(data_blue['inst_vega_mag1']-data_blue['inst_vega_mag2'], data_blue['inst_vega_mag2'], 'b.', markersize = 2)
        blim = info['faint_lim'][i] + 2
        ax4.set_ylim(blim, blim-6)
        ax4.set_xlim(-1, 3.5)
        ax4.set_xlabel('F606W-F814W')
        ax4.set_ylabel('F814W')
        ax4.set_title('CMD raw')
        
        color, mag = data_clip['inst_vega_mag1']-data_clip['inst_vega_mag2'], data_clip['inst_vega_mag2']
        judge_red, judge_blue = blue_selection(color, mag, info['red_lim'][i], info['faint_lim'][i])
        data_red, data_blue = data_clip[judge_red], data_clip[judge_blue]
        ax5 = fig.add_subplot(2, 18, (25, 30))
        ax5.plot(data_red['x'], data_red['y'], 'r.', markersize = 2)
        ax5.plot(data_blue['x'], data_blue['y'], 'b.', markersize = 2)
        ax5.set_xlabel('X')
        ax5.set_ylabel('Y')
        ax5.set_title('Spatial Distribution')
        
        ax6 = fig.add_subplot(2, 18, (33, 36))
        ax6.plot(data_clip['inst_vega_mag1']-data_clip['inst_vega_mag2'], data_clip['inst_vega_mag2'], 'k.', markersize = 2)
        ax6.set_xlim(-1, 3.5)
        ax6.set_xlabel('F606W-F814W')
        ax6.set_ylim(blim, blim-6)
        ax6.set_ylabel('F814W')
        ax6.set_title('CMD clipped')
        plt.savefig(os.path.join(os.path.split(__file__)[0], 'clipping_map/{:s}_{:s}.png'.format(info['field'][i], spclip)))
        plt.close()
        
spclip_list = ['5p', '10p', '20p']
for i1 in range(len(spclip_list)):
    spclip = spclip_list[i1]
    print(spclip)
    clipping(spclip)
'''
# determine color band
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

spclip_list = ['5p', '10p', '20p']
width_list = [0.75, 1, 1.5, 2]
slope_range_f606, slope_range_f555 = [-7, -4.95], [-3, -0.95]
info = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'sn_analysis_new.csv'))
for i in range(len(spclip_list)):
    clipname = spclip_list[i]
    for j in range(len(width_list)):
        width = width_list[j]
        print(clipname, width)
        slope, inter = np.zeros(len(info)), np.zeros(len(info))
        for k in range(len(info)):
            data = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'clipped_csv/{:s}_{:s}.csv'.format(info['field'][k], clipname)))
            color, mag = data['inst_vega_mag1']-data['inst_vega_mag2'], data['inst_vega_mag2']
            if (info['filter'][k] == 'F555W'):
                slope[k], inter[k] = color_band(color, mag, slope_range_f555, width)
            else:
                slope[k], inter[k] = color_band(color, mag, slope_range_f606, width)
        info['slope_ac'], info['inter_ac'] = slope, inter
        
        if width != 0.75:
            fname = 'clip={:s},width={:.1f}'.format(clipname, width)
        else:
            fname = 'clip={:s},width={:.2f}'.format(clipname, width)
        info.to_csv(os.path.join(os.path.split(__file__)[0], 'info/sn_info_{:s}.csv'.format(fname)), index=False)
