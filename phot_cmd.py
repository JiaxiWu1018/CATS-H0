"""
Created on Sat Feb 04 18:03 2023
@author: Jiaxi Wu
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

info = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'sn_analysis.csv'))
for i in range(len(info)):
    data = pd.read_csv(os.path.join(os.path.split(__file__)[0], 'phot_csv/{:s}.csv'.format(info['field'][i])))
    fig, ax = plt.subplots(1, 1, figsize=(5, 8))
    color, mag = data['inst_vega_mag1']-data['inst_vega_mag2'], data['inst_vega_mag2']
    ax.plot(color, mag, 'k.', ms=2)
    ax.set_xlim(-1, 3.5)
    ax.set_ylim(30, 22)
    ax.set_title(info['field'][i])
    plt.savefig(os.path.join(os.path.split(__file__)[0], 'phot_cmd/{:s}.png'.format(info['field'][i])), bbox_inches='tight')
    plt.close()