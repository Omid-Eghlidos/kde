#!/usr/bin/env python3
"""
This script processes structural distribution files generated in the current
directory and plots them.
It supports Radial Distribution Function (RDF), Bond Distribution Function (BDF),
Angle Distribution Function (ADF), Torsion Distribution Function (TDF), and
Improper Distribution Function (IDF).

For each distribution type, it reads the corresponding text files, applies a
Savitzky-Golay filter for smoothing,
converts angles to radians for angular distributions, and saves the plots in a
specified folder (default: 'plots').

Dependencies:
- numpy
- scipy
- matplotlib

The script is designed to be run in a directory containing the structural
distribution text files.

Usage:
    python3 plot_distributions.py

The generated plots will be saved in the 'plots' folder.
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
    """
    Main function to go through different structural distribution files in the
    current directory, generate plots for each, and save them in a specified folder.
    """
    path = make_folder('plots')
    dfs = ['rdf', 'bdf', 'adf', 'tdf', 'idf']
    for df in dfs:
        plot_distribution(path, df)


def plot_distribution(path, df):
    """
    Reads data from the specified distribution files, applies filtering,
    plots the distributions, and saves them as PNG images.

    Parameters:
    - path: str, the directory where plots will be saved.
    - df: str, the type of distribution (e.g., 'rdf', 'bdf').

    Returns:
    - None
    """
    files = glob(f'*{df}*.txt')
    names = dict(rdf='Pair', bdf='Bond', adf='Angle', tdf='Dihedral', idf='Improper')
    for f in files:
        Type = f.strip(df).strip('.txt').split('_')[-1]
        data = np.genfromtxt(f)
        data[:,1] = savgol_filter(data[:,1], 11, 2)
        if df in ['adf', 'tdf', 'idf']:
            data[:,0] *= np.pi / 180.0
            data[:,1] *= 180.0 / np.pi

        fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.8))
        ax.plot(data[:, 0], data[:, 1], label=f'{names[df]} {Type}')
        adjust_plot_settings(df, ax)
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(f'{path}/{df}_{Type}.png', dpi=300, bbox_inches='tight')
        plt.close()


def adjust_plot_settings(df, ax):
    """
    Configures axis labels, limits, and ticks based on the distribution type.

    Parameters:
    - df: str, the type of distribution (e.g., 'rdf', 'bdf').
    - ax: matplotlib.axes.Axes, the axis object to configure.

    Returns:
    - None
    """
    if df == 'rdf':
        label = dict(x='$r\;(\AA)$', y='$g(r)$')
        ax.set_xlim(0.0, 15.0)
    elif df == 'bdf':
        label = dict(x='$l\;(\AA)$', y='$P(l)\;(\AA^{-1})$')
        ax.set_xlim(0.0, 5.0)
    elif df == 'adf':
        label = dict(x='$\\theta$ (rad)', y='$P(\\theta)$ (rad$^{-1}$)')
        xmin, xmax = 0.0, np.pi
        ax.set_xlim(xmin, xmax)
        xticks = np.linspace(xmin, xmax, 7)
        ax.set_xticks(xticks)
        ax.set_xticklabels(['$0$',  '$\pi/6$',  '$\pi/3$', '$\pi/2$',
                                   '$2\pi/3$', '$5\pi/6$',   '$\pi$'])
    elif df == 'tdf':
        label = dict(x='$\phi$ (rad)', y='$P(\\phi)$ (rad$^{-1}$)')
        xmin, xmax = -np.pi, np.pi
        ax.set_xlim(xmin, xmax)
        xticks = np.linspace(xmin, xmax, 7)
        ax.set_xticks(xticks)
        ax.set_xticklabels([ '$-\pi$', '$-2\pi/3$', '$-\pi/3$', '$0$',
                            '$\pi/3$',  '$2\pi/3$',   '$\pi$'])
    elif df == 'idf':
        label = dict(x='$\psi$ (rad)', y='$P(\\psi)$ (rad$^{-1})$')
        xmin, xmax = 0.0, np.pi
        ax.set_xlim(xmin, xmax)
        xticks = np.linspace(xmin, xmax, 7)
        ax.set_xticks(xticks)
        ax.set_xticklabels(['$0$',  '$\pi/6$',  '$\pi/3$', '$\pi/2$',
                                   '$2\pi/3$', '$5\pi/6$',   '$\pi$'])

    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(label['x'])
    ax.set_ylabel(label['y'])


def make_folder(name):
    """
    Creates a folder with the specified name in the current working directory
    if it does not already exist.

    Parameters:
    - name: str, the name of the folder to create.

    Returns:
    - str: The full path to the created (or existing) folder.
    """
    path = os.path.join(os.getcwd(), name)
    if not os.path.exists(path):
        os.mkdir(path)
    return path


if __name__ == '__main__':
    main()

