"""
Created on Sun Dec 18 21:19 2022
@author: Jiaxi Wu
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import linalg
import os

# second part - TRGB detection
def selection(xx, yy, band, width):
    m, b = band[0], band[1]
    xx1 = (yy - b) / m
    xx2 = (yy - b) / m + width
    judge = (xx > xx1) & (xx < xx2)
    xx, yy = xx[judge], yy[judge]
    return xx, yy

def gloess(mag, tau, bin_width = 0.01):
    fig_gloess, ax_gloess = plt.subplots(figsize = (10, 10))
    hist, bins, cont = ax_gloess.hist(mag, bins = np.arange(min(mag), max(mag) + bin_width, bin_width), color = 'black')
    plt.close()
    bin_centers = []
    for j in range(len(bins) - 1):
        bin_centers.append((bins[j] + bins[j + 1]) / 2)
    yest = np.zeros(len(hist))
    w = np.array([np.exp(- (bin_centers - bin_centers[i])**2/(2 * tau**2)) for i in range(len(hist))])
    for i in range(len(hist)):
        weights = w[i, :]
        b = np.array([np.sum(weights * hist), np.sum(weights * hist * bin_centers)])
        A = np.array([[np.sum(weights), np.sum(weights * bin_centers)],
                    [np.sum(weights * bin_centers), np.sum(weights * bin_centers * bin_centers)]])
        theta = linalg.solve(A, b)
        yest[i] = theta[0] + theta[1] * bin_centers[i]
    return yest, np.array(bin_centers[2: -2])

def poisson_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append(hist[ii - 1] * -1 * 1/np.sqrt(hist[ii - 1]) + hist[ii] * 0 + hist[ii + 1] * 1* 1/np.sqrt(hist[ii + 1]))
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def hatt_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append((hist[ii + 1] - hist[ii - 1]) / np.sqrt(hist[ii + 1] + hist[ii - 1]))
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def simple_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append(hist[ii-1] * -1 + hist[ii+1])
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def prod_edr(mag, tau, weighting):
    hist, binval = gloess(mag, tau)
    if (weighting == 'poisson'):
        lumfunc, edres = poisson_sobel(hist)
    elif (weighting == 'hatt'):
        lumfunc, edres = hatt_sobel(hist)
    elif (weighting == 'simple'):
        lumfunc, edres = simple_sobel(hist)
    else:
        print('No such weighting type')
        return
    return mag, binval, edres

def qualify(trgb, mag, ratio_min, nbt_min):
    RGB_pop = sum((mag < trgb + 0.5) & (mag > trgb))
    AGB_pop = max(sum((mag < trgb) & (mag > trgb - 0.5)), 1)
    nbt = sum((mag > trgb) & (mag < trgb+1))
    R = RGB_pop / AGB_pop
    judge = (R >= ratio_min) & (nbt >= nbt_min)
    return judge, R, nbt

def findpeak(binval, func, mag, ratio_min, nbt_min, peak_min):
    pos, peak, fwhm, ratio, nbt, peakq = [], [], [], [], [], []
    for i in range(1, len(func) - 1):
        if (func[i] > func[i - 1]) & (func[i] > func[i + 1]):
            left_edge, right_edge = i, i
            while (func[left_edge] > func[left_edge - 1]) & (func[left_edge] > func[i] / 2) & (left_edge > 0):
                left_edge -= 1
            while (func[right_edge] > func[right_edge + 1]) & (func[right_edge] > func[i] / 2) & (right_edge < len(func) - 2):
                right_edge += 1
            
            test, R, NBT = qualify(binval[i], mag, ratio_min, nbt_min)
            pos.append(i)
            ratio.append(R)
            nbt.append(NBT)
            peak.append(func[i])
            fwhm.append((right_edge - left_edge) / 100)
            if test == True:
                peakq.append(func[i])

    pos_good, ratio_good, nbt_good, fwhm_good = [], [], [], []
    EDR, percent = [], []
    if len(peakq) > 0:
        maxv = max(peakq)
        for i in range(len(pos)):
            if peak[i] >= peak_min * maxv:
                pos_good.append(pos[i])
                ratio_good.append(ratio[i])
                nbt_good.append(nbt[i])
                fwhm_good.append(fwhm[i])
                EDR.append(peak[i])
                percent.append(peak[i]/maxv)
    trgb = binval[pos_good]
    return trgb, np.array(ratio_good), np.array(nbt_good), np.array(fwhm_good), np.array(EDR), np.array(percent)

def cal_err(R, N):
    sigma = np.sqrt((2*np.exp(1.5*(3-R))/(np.exp(1.5*(3-R))+1)*(1/(N-100))**0.1)**2+0.04**2)
    return sigma

def detection(spclip, width, smoothing, peak_min, int_ext, weighting):
    ratio_min, nbt_min = 3, 0
    ratio_cut, nbt_cut = 1.5, 101
    band = 'yes'
    if width != 0.75:
        infoname = 'clip={:s},width={:.1f}'.format(spclip, width)
    else:
        infoname = 'clip={:s},width={:.2f}'.format(spclip, width)
    fname = '{:s},smooth={:.2f},min_th={:.1f},int_ext={:s}'.format(infoname, smoothing, peak_min, int_ext)
    info = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'info/sn_info_{:s}.csv'.format(infoname)))

    agg_galaxy, agg_galra, agg_galdec = np.array([]), np.array([]), np.array([])
    agg_field, agg_fdra, agg_fddec = np.array([]), np.array([]), np.array([])
    agg_trgb, agg_ratio, agg_raerr, agg_ratio2 = np.array([]), np.array([]), np.array([]), np.array([])
    agg_nbt, agg_EDR, agg_percent = np.array([]), np.array([]), np.array([])
    agg_fwhm, agg_slope, agg_color = np.array([]), np.array([]), np.array([])
    agg_a606, agg_a814, agg_ext = np.array([]), np.array([]), np.array([])
    agg_inc, agg_pos, agg_d2n = np.array([]), np.array([]), np.array([])

    for i in range(len(info)):
        data = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'clipped_csv/{:s}_{:s}.csv'.format(info['field'][i], spclip)))
        color, mag = data['inst_vega_mag1']-data['inst_vega_mag2'], data['inst_vega_mag2']

        band = np.array([info['slope_ac'][i], info['inter_ac'][i]])
        color, mag = selection(color, mag, band, width)
        
        mag, binval, edres = prod_edr(mag, smoothing, weighting)
        trgb, ratio, nbt, fwhm, EDR, percent = findpeak(binval, edres, mag, ratio_min, nbt_min, peak_min)
        judge = (ratio >= ratio_cut) & (nbt >= nbt_cut)
        trgb, ratio, nbt, fwhm = trgb[judge], ratio[judge], nbt[judge], fwhm[judge]
        EDR, percent = EDR[judge], percent[judge]
        
        tipcol, tipsl, raerr, ra2 = np.array([]), np.array([]), np.array([]), np.array([])
        judge = [True]*len(trgb)
        for k in range(len(trgb)):
            judgetip = (mag < trgb[k]+0.1) & (mag > trgb[k]-0.1)
            judge1b = (mag < trgb[k]+1.1) & (mag > trgb[k]+0.9)
            if (np.sum(judgetip) == 0):
                judge[k]=False
                continue
            coltip = np.mean(color[judgetip])
            tipcol = np.append(tipcol, np.array([coltip]))
            if (np.sum(judge1b) == 0):
                tipsl = np.append(tipsl, np.array([-99]))
            else:
                col1b = np.mean(color[judge1b])
                tipsl = np.append(tipsl, np.array([1/(coltip-col1b)]))
            RGB = np.sum((mag>trgb[k])&(mag<trgb[k]+0.5))
            AGB = max(np.sum((mag<trgb[k])&(mag>trgb[k]-0.5)),1)
            tipra = RGB/AGB
            tipraerr = tipra * np.sqrt(1/RGB + 1/AGB)
            raerr = np.append(raerr, np.array([tipraerr]))
            RGB2 = np.sum((mag>trgb[k])&(mag<trgb[k]+1))
            AGB2 = max(np.sum((mag<trgb[k])&(mag>trgb[k]-1)),1)
            ra2 = np.append(ra2, np.array([RGB2/AGB2]))
        trgb, ratio, nbt, fwhm = trgb[judge], ratio[judge], nbt[judge], fwhm[judge]
        EDR, percent = EDR[judge], percent[judge]
        
        agg_trgb, agg_ratio, agg_nbt = np.append(agg_trgb, trgb), np.append(agg_ratio, ratio), np.append(agg_nbt, nbt)
        agg_raerr, agg_ratio2 = np.append(agg_raerr, raerr), np.append(agg_ratio2, ra2)
        agg_slope, agg_color = np.append(agg_slope, tipsl), np.append(agg_color, tipcol)
        agg_fwhm, agg_EDR, agg_percent = np.append(agg_fwhm, fwhm), np.append(agg_EDR, EDR), np.append(agg_percent, percent)
        agg_galaxy = np.append(agg_galaxy, np.array([info['galaxy'][i]]*len(trgb)))
        agg_galra = np.append(agg_galra, np.array([info['gal_ra'][i]]*len(trgb)))
        agg_galdec = np.append(agg_galdec, np.array([info['gal_dec'][i]]*len(trgb)))
        agg_field = np.append(agg_field, np.array([info['field'][i]]*len(trgb)))
        agg_fdra = np.append(agg_fdra, np.array([info['fd_ra'][i]]*len(trgb)))
        agg_fddec = np.append(agg_fddec, np.array([info['fd_dec'][i]]*len(trgb)))
        agg_ext = np.append(agg_ext, np.array([info['int_ext'][i]]*len(trgb)))
        agg_a606 = np.append(agg_a606, np.array([info['A555/606'][i]]*len(trgb)))
        agg_a814 = np.append(agg_a814, np.array([info['A814'][i]]*len(trgb)))
        agg_inc = np.append(agg_inc, np.array([info['inc'][i]]*len(trgb)))
        agg_pos = np.append(agg_pos, np.array([info['pos_ang'][i]]*len(trgb)))
        agg_d2n = np.append(agg_d2n, np.array([info['dist2nuc'][i]]*len(trgb)))

    data = pd.DataFrame()
    data['galaxy'], data['gal_ra'], data['gal_dec'] = agg_galaxy, agg_galra, agg_galdec
    data['field'], data['fd_ra'], data['fd_dec'] = agg_field, agg_fdra, agg_fddec
    data['TRGB'], data['TRGB_err'] = agg_trgb, cal_err(agg_ratio, agg_nbt)
    data['R'], data['R_err'], data['R_pm1'], data['N'] = agg_ratio, agg_raerr, agg_ratio2, agg_nbt
    data['EDR_value'], data['percentage'], data['FWHM'], data['slope'], data['tip_color'] = agg_EDR, agg_percent, agg_fwhm, agg_slope, agg_color
    data['A555/606'], data['A814'], data['int_ext'] = agg_a606, agg_a814, agg_ext
    data['inc'], data['pos_ang'], data['dist2nuc'] = agg_inc, agg_pos, agg_d2n

    if int_ext == 'no':
        data['int_ext'] = np.zeros(len(data))

    data = data.sort_values(by=['galaxy', 'field'], ascending=[True, True])
    data = data.reset_index(drop=True)
    data.to_csv(os.path.join(os.path.split(__file__)[0], 'detection/sn_detection_{:s}.csv'.format(fname)), index = False)

# set parameter
spclip_list = ['20p'] #['5p', '10p', '20p']
width_list = [0.75, 1.0, 1.5, 2.0]
smoothing_list = [0.07, 0.10, 0.15]
peak_min_list = [0.5, 0.6, 0.8]
int_ext_list = ['yes', 'no']
#weighting_list = ['simple', 'hatt']

for i1 in range(len(spclip_list)):
    spclip = spclip_list[i1]
    for i2 in range(len(width_list)):
        width = width_list[i2]
        for i3 in range(len(smoothing_list)):
            smoothing = smoothing_list[i3]
            for i4 in range(len(peak_min_list)):
                peak_min = peak_min_list[i4]
                for i5 in range(len(int_ext_list)):
                    int_ext = int_ext_list[i5]
                    print(spclip, width, smoothing, peak_min, int_ext)
                    detection(spclip, width, smoothing, peak_min, int_ext, weighting='hatt')