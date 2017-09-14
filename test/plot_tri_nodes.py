"""
Plot the SBP tri nodes
"""
dump=False # set to True to write a png file
facecub=False #True # include the face cubature points?
degree = 4

import matplotlib
if dump:
    #matplotlib.use('PS')
    matplotlib.use('Agg') # for PNG

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math

# set some formating parameters
axis_fs = 10 # axis title font size
axis_lw = 1.0 # line width used for axis box, legend, and major ticks
label_fs = 8 # axis labels' font size

# set figure size in inches, and crete a single set of axes
fig = plt.figure(figsize=(2,np.sqrt(3)), facecolor='w', dpi=300)
ax = fig.add_subplot(111)

data = open("./nodes.dat")
#data = open("./p" + str(degree) + "_bndry.dat", 'r')
#data = open("./p" + str(degree) + "_intr.dat", 'r')
arrays = [np.array(map(float, line.split())) for line in data]
x = arrays[0]
y = arrays[1]
#xf = arrays[2]
#yf = arrays[3]
#x, y, xf, yf = np.loadtxt(data)

#tri, = ax.plot(np.append(x[0:3],x[0]), np.append(y[:3],y[0]), 'k-', lw=1)
tri, = ax.plot([0.0, 1.0, 0.5, 0.0], [0.0, 0.0, np.sqrt(3)*0.5, 0.0], '-k', lw=1)
#nodes, = ax.plot(x, y, 'ko', lw=1, mfc='w', ms=6, mec='k', mew=1.0)
plt.gca().set_aspect('equal')

label = False
if label:
    idx = 1
    props = dict(boxstyle='round', facecolor='white')
    for coord in zip(x,y):
        ax.text(coord[0], coord[1],  idx, \
                size=label_fs, ha='center', \
                horizontalalignment='center', verticalalignment='center', bbox=props)
        idx += 1
else:
    if facecub:
        facenodes, = ax.plot(xf, yf, 'ks', lw=1, ms=5, mec='k', mew=1.0)
    nodes, = ax.plot(x, y, 'ko', lw=1, mfc='w', ms=5, mec='k', mew=1.0)

if 1 == 0:
    # Tweak the appeareance of the axes
    #ax.axis([0., 1., 0., 1.])  # axes ranges
    ax.axis([-0.1, 1.1, -0.1, np.sqrt(3)*0.5+0.1])  # axes ranges
    ax.set_position([0.15, 0.15, 0.8, 0.8]) # position relative to figure edges
    #ax.set_xlabel(r"$\xi$", fontsize=axis_fs, weight='bold', labelpad=0)
    ax.set_xlabel(r"$x$", fontsize=axis_fs, weight='bold', labelpad=0)
    ax.xaxis.set_label_coords(0.5, -0.12)
    #ax.set_ylabel(r"$\eta$", fontsize=axis_fs, weight='bold', rotation=0)
    ax.set_ylabel(r"$y$", fontsize=axis_fs, weight='bold', rotation=0)
    ax.yaxis.set_label_coords(-0.17, 0.5)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(axis_lw)

    # ticks on bottom and left only
    ax.xaxis.tick_bottom() # use ticks on bottom only
    ax.yaxis.tick_left()
    for line in ax.xaxis.get_ticklines():
        line.set_markersize(-6) # length of the tick
        line.set_markeredgewidth(axis_lw) # thickness of the tick
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(-6) # length of the tick
        line.set_markeredgewidth(axis_lw) # thickness of the tick
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(label_fs)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(label_fs)
    for tick in ax.get_xaxis().get_major_ticks():
        tick.set_pad(8.)
        tick.label1 = tick._get_text1()
    for tick in ax.get_yaxis().get_major_ticks():
        tick.set_pad(8.)
        tick.label1 = tick._get_text1()

    # We change the fontsize of minor ticks label
    ax.tick_params(axis='both', which='major', labelsize=label_fs)
    ax.tick_params(axis='both', which='minor', labelsize=label_fs)
    #ax.axhline(linewidth=axis_lw)
    #ax.axvline(linewidth=axis_lw)
else:
    #ax.axis([-0.1, 1.1, -0.1, np.sqrt(3)*0.5+0.1])  # axes ranges
    #ax.set_position([0.01, 0.01, 0.99, 0.99]) # position relative to figure edges
    ax.set_position([0.0, 0.0, 1., 1.]) # position relative to figure edges
    plt.axis("off")
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

if dump:
    #fig.savefig('temp.eps', facecolor=fig.get_facecolor(), dpi=300, edgecolor='none')
    fig.savefig('temp.png', facecolor=fig.get_facecolor(), dpi=300, edgecolor='none')
else:
    plt.show()
