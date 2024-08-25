#!/usr/bin/env python3
from glob import glob
import numpy
from matplotlib import pyplot
import scipy.stats

# Set default plotting parameters
pyplot.rcParams['font.family'] = 'serif'
pyplot.rcParams['font.serif'] = ['Times New Roman']
pyplot.rcParams['font.size'] = 9
pyplot.rcParams['mathtext.fontset'] = 'stix'

def plot_rdfs(tag):
    # Read data
    data = []
    for f in glob('*rdf*.txt'):
        data.append(numpy.genfromtxt(f))

    # Compute statistics (mean and confidence intervals.)
    rdf = numpy.mean(data, axis=0)

    # Write the average into file
    fid = open("rdf-{}".format(tag), 'w')
    for i in range(len(rdf)):
        fid.write("{:6.2f}\t{:6.3f}\n".format(rdf[i,0], rdf[i,1]))

    # T statistic for 90% confidence interval.
    t_value = scipy.stats.t.ppf(1.0-0.05, len(data)-1)
    ci = numpy.std(data, axis=0)[:,1] * t_value / numpy.sqrt(len(data))

    # Plot
    pyplot.clf()
    pyplot.figure(figsize=(4,3))
    ax = pyplot.subplot(111)
    #ax.fill_between(rdf[:,0], rdf[:,1]-ci, rdf[:,1]+ci, alpha=0.2, lw=0)
    ax.plot(rdf[:,0], rdf[:,1], label='RDF {}'.format(tag))
    ax.set_xlim(0.0, 14.0)
    #ax.set_ylim(0.0, 2.5)
    ax.set_xlabel('$r$ ($\AA$)')
    ax.set_ylabel('$g(r)$')
    pyplot.legend(loc='best')
    pyplot.tight_layout()
    pyplot.savefig('rdf-{}.png'.format(tag), dpi=300)
    pyplot.close()


def plot_bdfs(tag):
    # Reda data
    data = []
    for f in glob('*bdf*.txt'):
        data.append(numpy.genfromtxt(f))

    # Compute statistics (mean and confidence intervals.)
    bdf = numpy.mean(data, axis=0)
    
    # Write the average into file
    fid = open("bdf-{}".format(tag), 'w')
    for i in range(len(bdf)):
        fid.write("{:6.2f}\t{:6.3f}\n".format(bdf[i,0], bdf[i,1]))

    # T statistic for 90% confidence interval.
    t_value = scipy.stats.t.ppf(1.0-0.05, len(data)-1)
    ci = numpy.std(data, axis=0)[:,1] * t_value / numpy.sqrt(len(data))

    # Plot
    pyplot.clf()
    pyplot.figure(figsize=(4,3))
    ax = pyplot.subplot(111)
    #ax.fill_between(bdf[:,0], bdf[:,1]-ci, bdf[:,1]+ci, alpha=0.2, lw=0)
    ax.plot(bdf[:,0], bdf[:,1], label='BDF {}'.format(tag))
    ax.set_xlim(0.0, 4.0)
    #ax.set_ylim(0, 13.0)
    ax.set_xlabel('$l$ ($\AA$)')
    ax.set_ylabel('$P(l)$')
    pyplot.legend(loc='best')
    pyplot.tight_layout()
    pyplot.savefig('bdf-{}.png'.format(tag), dpi=300)
    pyplot.close()


def plot_adfs(tag):
    # Read data
    data = []
    for f in glob('*adf*.txt'):
        data.append(numpy.genfromtxt(f))

    # Compute statistics (mean and confidence intervals.)
    adf = numpy.mean(data, axis=0)
    
    # Write the average into file
    fid = open("adf-{}".format(tag), 'w')
    for i in range(len(adf)):
        fid.write("{:6.2f}\t{:6.3f}\n".format(adf[i,0], adf[i,1]))

    # T statistic for 90% confidence interval.
    t_value = scipy.stats.t.ppf(1.0-0.05, len(data)-1)
    ci = numpy.std(data, axis=0)[:,1] * t_value / numpy.sqrt(len(data))

    # Plot
    pyplot.clf()
    pyplot.figure(figsize=(4,3))
    ax = pyplot.subplot(111)
    #ax.fill_between(adf[:,0], adf[:,1]-ci, adf[:,1]+ci, alpha=0.2, lw=0)
    ax.plot(adf[:,0], adf[:,1], label='ADF {}'.format(tag))
    ax.set_xlim(0.0, 180.0)
    #ax.set_ylim(0.0, 0.10)
    ax.set_xlabel('$\\theta^\\circ$')
    ax.set_ylabel('$P$($\\theta$)')
    pyplot.legend(loc='best')
    pyplot.tight_layout()
    pyplot.savefig('adf-{}.png'.format(tag), dpi=300)
    pyplot.close()


def plot_tdfs(tag):
    # Read data
    data = []
    for f in glob('*tdf*.txt'):
        data.append(numpy.genfromtxt(f))

    # Compute statistics (mean and confidence intervals.)
    tdf = numpy.mean(data, axis=0)
    
    # Write the average into file
    fid = open("tdf-{}".format(tag), 'w')
    for i in range(len(tdf)):
        fid.write("{:6.2f}\t{:6.3f}\n".format(tdf[i,0], tdf[i,1]))

    # T statistic for 90% confidence interval.
    t_value = scipy.stats.t.ppf(1.0-0.05, len(data)-1)
    ci = numpy.std(data, axis=0)[:,1] * t_value / numpy.sqrt(len(data))

    # Plot
    pyplot.clf()
    pyplot.figure(figsize=(4,3))
    ax = pyplot.subplot(111)
    #ax.fill_between(tdf[:,0], tdf[:,1]-ci, tdf[:,1]+ci, alpha=0.2, lw=0)
    ax.plot(tdf[:,0], tdf[:,1], label='TDF {}'.format(tag))
    ax.set_xlim(-180.0, 180.0)
    #ax.set_ylim(0.0, 0.012)
    ax.set_xlabel('$\\phi^\\circ$')
    ax.set_ylabel('$P$($\\phi$)')
    pyplot.legend(loc='best')
    pyplot.tight_layout()
    pyplot.savefig('tdf-{}.png'.format(tag), dpi=300)
    pyplot.close()


def define_pair(tag):
    ''' Change the numeric tag to atom types '''
    pair = ""

    for i, t in enumerate(tag):
        if int(t) == 1:
            pair += "A"
        elif int(t) == 2:
            pair += "B"
        elif int(t) == 3:
            pair += "C"
        elif int(t) == 4:
            pair += "D"
        elif int(t) == 5:
            pair += "E"
        elif int(t) == 6:
            pair += "F"

        if i < len(tag) - 1:
            pair += "-"

    return pair


if __name__ == '__main__':
    plot_rdfs('AA')

    plot_bdfs('AA')

    plot_adfs('AAA')
   
    plot_tdfs('AAAA')
   

