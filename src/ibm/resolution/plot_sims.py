#!/usr/bin/env python3

import sys, re, os.path
import pandas as pd
import matplotlib

import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

f = open(filename);
fl = f.readlines();
f.close()


parline = -1

for idx, line in enumerate(fl):
    if re.match("^patch.*",line) != None:
        parline = idx - 1;
        break;

# read in the csv file
if parline > 0:
    histdat = pd.read_csv(filename, nrows=parline-3, sep=";")
else:
    histdat = pd.read_csv(filename, sep=";")

# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 3

# add first subplot
ax1 = plt.subplot(num_rows,1,1)
ax1.plot(
        histdat["generation"],histdat["s1"],'b',
        histdat["generation"],histdat["s2"],'r',
        histdat["generation"],histdat["qS"],'y',
        histdat["generation"],histdat["qNS"],'k'
        )
ax1.set_ylim([0,1])
ax1.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
ax1.set_ylabel(r'Phenotype')
ax1.legend((r'$s_{1}$',r'$s_{2}$',r'$q_{S}$',r'$q_{NS}$'))

# add second subplot
ax2 = plt.subplot(num_rows,1,2)
ax2.plot(histdat["generation"],histdat["prob1"],'b',
            histdat["generation"],histdat["prob2"],'r')
ax2.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
ax2.set_ylabel(r'Prop $z_{1}$')
ax2.set_ylim([0,1])
ax2.legend((r'$p_{1}$',r'$p_{2}$'))

# add third subplot
plt.subplot(num_rows,1,3)
plt.plot(
        histdat["generation"],histdat["vars1"],'b',
        histdat["generation"],histdat["vars2"],'r',
        histdat["generation"],histdat["varqNS"],'g',
        histdat["generation"],histdat["varqS"],'black',
        linewidth=1)
#plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Variance')
plt.legend((r'$\sigma_{s_{1}}^{2}$',r'$\sigma_{s_{2}}^{2}$',r'$\sigma_{q_{NS}}^{2}$',r'$\sigma_{q_{S}}^{2}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
