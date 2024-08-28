#!/usr/bin/env python3
"""
Goal
----
This code goes over all the generated structural distributions text file in the
current folder and plot and save them in a folder named "plots".
"""

import os
import numpy as np
from glob import glob
from scipy.signal import savgol_filter
from matplotlib import pyplot as plt

plt.rc('font', family = 'Times New Roman')
plt.rc('font', size = 8)
plt.rc('mathtext', fontset = 'stix')


def main():
    ''' Go through different distributions and plot them if their file exists. '''
    path = make_folder('plots')
    dfs = ['rdf', 'bdf', 'adf', 'tdf', 'idf']
    for df in dfs:
        plot_distribution(path, df)


def plot_distribution(path, df):
    ''' Read data from all the distribution files, plot, and store them in a
        separate folder. '''
    files = glob(f'*{df}*.txt')
    for f in files:
        Type = f.strip(df).strip('.txt').split('_')[-1]
        data = np.genfromtxt(f)
        data[:,1] = savgol_filter(data[:,1], 11, 2)
        if df in ['adf', 'tdf', 'idf']:
            data[:,1] *= 180.0 / np.pi

        fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))
        label, lim = adjust_plot_settings(df)
        ax.plot(data[:, 0], data[:, 1], label=f'{label["plot"]} {Type}')
        ax.set_xlim(lim['x'][0], lim['x'][1])
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel(label['x'])
        ax.set_ylabel(label['y'])
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(f'{path}/{df}_{Type}.png', dpi=300, bbox_inches='tight')
        plt.close()


def adjust_plot_settings(df):
    ''' Return axis and label settings for the given structural distributions. '''
    if df == 'rdf':
        label = dict(x='$r\;(\AA)$', y='$g(r)$', plot='Pair')
        lim = dict(x=[0.0, 15.0], y=[0.0, 15.0])
    elif df == 'bdf':
        label = dict(x='$l\;(\AA)$', y='$P(l)\;(\AA^{-1})$', plot='Bond')
        lim = dict(x=[0.0, 5.0], y=[0.0, 1.0])
    elif df == 'adf':
        label = dict(x='$\\theta$ (rad)', y='$P(\\theta)$ (rad$^{-1}$)', plot='Angle')
        lim = dict(x=[0.0, 180.0], y=[0.0, 0.1])
    elif df == 'tdf':
        label = dict(x='$\phi$ (rad)', y='$P(\\phi)$ (rad$^{-1}$)', plot='Dihedral')
        lim = dict(x=[-180.0, 180.0], y=[0.0, 0.1])
    elif df == 'idf':
        label = dict(x='$\psi$ (rad)', y='$P(\\psi)$ (rad$^{-1})$', plot='Improper')
        lim = dict(x=[0.0, 180.0], y=[0.0, 0.1])
    return label, lim



def make_folder(dir_name):
    ''' Make folder in the current working directory, if it does not exist. '''
    path = os.path.join(os.getcwd(), dir_name)
    try:
        os.mkdir(path)
        exist = False
    except OSError:
        exist = True
    return path


if __name__ == '__main__':
    main()

