"""
Plot the SBP nodes
"""
dump=False # set to True to write a png file
import matplotlib
if dump:
    matplotlib.use('Agg')

from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import PolyCollection
#from matplotlib.colors import colorConverter
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math

# set some formating parameters
axis_fs = 12 # axis title font size
axis_lw = 1.0 # line width used for axis box, legend, and major ticks
label_fs = 8 # axis labels' font size

# set figure size in inches, and crete a single set of axes
fig = plt.figure(figsize=(2.0,2.0), facecolor='w', dpi=300)
ax = fig.add_subplot(111, projection='3d', aspect='equal')

#data = open("./p" + str(degree) + "_intr.dat", 'r')
data = open("./nodes.dat", 'r')
#data = open("./nodes_variant7.dat", 'r')
#data = open("./p2_nodes_variant2.dat", 'r')
#data = open("./p4_nodes_variant6.dat", 'r')
arrays = [np.array(map(float, line.split())) for line in data]
x = arrays[0]
y = arrays[1]
z = arrays[2]
#xf = arrays[3]
#yf = arrays[4]
#zf = arrays[5]

#degree = 1

x_vtx = np.array([0.0,1.0,0.5,0.5])
y_vtx = np.array([0.0,0.0,np.sqrt(3)/2,np.sqrt(3)/6])
z_vtx = np.array([0.0,0.0,0.0,np.sqrt(2.0/3.0)])

ax.plot(np.append(x_vtx[0:4],[x_vtx[0],x_vtx[2],x_vtx[3],x_vtx[1]]),\
        np.append(y_vtx[0:4],[y_vtx[0],y_vtx[2],y_vtx[3],y_vtx[1]]), \
        np.append(z_vtx[0:4],[z_vtx[0],z_vtx[2],z_vtx[3],z_vtx[1]]), \
        '-k', lw=1, mfc='w', ms=4, mec='k', mew=1.0)
#ax.plot(x, y, z, 'ko', lw=1, mfc='w', ms=4, mec='k', mew=1.0)  #marker='o', color='w', s=4) #, depthshade=False)

ax.scatter(x, y, z, color='k', s=40, linewidths=0.0, alpha=1.0, depthshade=True)

#ax.scatter(xf, yf, zf, color='r', s=20, linewidths=0.0, alpha=1.0)

# verts = []
# verts.append(list(zip(x[0:3],y[0:3]))) #,np.array(z[0:3])))
# print verts
# cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.6)
# poly = PolyCollection(verts) #, facecolors = [cc('k')]) #, cc('g'), cc('b'),cc('y')])
# poly.set_alpha(0.7)
# ax.add_collection3d(poly) #, zs=zs, zdir='z')

if 1 == 0:
    ax.plot_trisurf(x[0:3], y[0:3], z[0:3], cmap=cm.gray, shade=True, linewidth=1, alpha=0.3)
    ax.plot_trisurf(x[1:4], y[1:4], z[1:4], cmap=cm.gray, shade=True, linewidth=1, alpha=0.3)
    ax.plot_trisurf(np.array([x[0],x[2],x[3]]), \
                np.array([y[0],y[2],y[3]]), \
                np.array([z[0],z[2],z[3]]), \
                cmap=cm.gray, shade=True, linewidth=1, alpha=0.3)
    ax.plot_trisurf(np.array([x[0],x[1],x[3]]), \
                np.array([y[0],y[1],y[3]]), \
                np.array([z[0],z[1],z[3]]), \
                cmap=cm.gray, shade=True, linewidth=1, alpha=0.1)

ax.set_position([0.15, 0.1, 1., 1.0]) # position relative to figure edges
plt.axis("off")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

ax.autoscale_view(tight=True, scaley=True, scalez=True) #, scalex=True, scaley=True, scalez=True)
#ax.zoom(1.0)
print ax.can_pan()
ax.view_init(elev=20, azim=-50)
ax.dist=6 # distance from camera to origin?


for o in fig.findobj():
    o.set_clip_on(False)

if dump:
    fig.savefig('temp.pdf', facecolor=fig.get_facecolor(), dpi=300, edgecolor='none')
    #fig.savefig('temp.png', facecolor=fig.get_facecolor(), dpi=300, edgecolor='none')
else:
    plt.show()
