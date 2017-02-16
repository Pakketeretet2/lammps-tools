import gen_morphexpt as gm
import numpy as np
import scipy
import sys
import os, glob
import re
import skeleton_tools as skt
import cPickle as pickle
from mayavi import mlab
import networkx as nx

def plot_3d(coords, figidx = None, color = None, colormap = None, vmin = None, vmax = None, scale_factor = 1., resolution = 15, connection = False, adjacent = False, bondcolor = (1, 1, 1),particle_diameter=1.0, neighborcutoff=0.2,scalemode='none'):
    '''
    3D visualization of particles. with option of showing Delaunay triangulations.
    '''
    if figidx == None:
        f = mlab.figure(size = (800, 838))
    else:
        f = mlab.figure(figidx)
    if color is None:
        color = tuple(np.random.random(3))
        points=mlab.points3d(coords[:,0], coords[:,1], coords[:,2], color = color, resolution=resolution, 
                             scale_factor=scale_factor, scale_mode='none')
    else:
        if colormap is None: colormap = 'jet'
        if vmin is None: vmin = color.min()
        if vmax is None: vmax = color.max()
    points=mlab.points3d(coords[:,0], coords[:,1], coords[:,2], color, colormap=colormap, resolution= resolution, 
                             scale_factor=scale_factor, scale_mode=scalemode, vmin = vmin, vmax = vmax )
    if connection is True:
        connections = np.array(delaunay_graph(coords).edges())
        points.mlab_source.dataset.lines = connections
        points.mlab_source.update()
        mlab.pipeline.surface(points, color=bondcolor, representation='wireframe')
    if adjacent is True:
        adjm = adjacentmatrix(coords,neighborcutoff=neighborcutoff,particle_diameter=particle_diameter)
        g = nx.Graph(adjm)
        connections = np.array(g.edges())
        points.mlab_source.dataset.lines = connections
        points.mlab_source.update()
        mlab.pipeline.surface(points, color=bondcolor, representation='wireframe')    
    return points



ptfile = sys.argv[1]
print "file to parse:", ptfile

# load the original point data, from a file with three columns depicting x,y,z coordinates of points# points are assumed to have diameter of 1
mlab.figure(2)
pts = np.loadtxt(ptfile)
plot_3d(pts,resolution=12,figidx=2)



# convert into a MorphExpt object, which merges the points into a continous solid
# and performs various morphological operations that can be used to extract
# the average width etc.
me = gm.gen_morphexpt(pts)

# this plot shows the gridded data generated by the morphology.
# each point is coloured by its distance from the closest edge of the domain.
mlab.figure(1)
plot_3d(me.domain_pts,color=me.sphere_distance(),resolution=3,scale_factor = 1./8,figidx=1)

# plot the approximation to the morphological skeleton of the domain

# Each skeleton is a list of points that traces out the path in 3D, and a list of distances
# of those points from the domain edge. The example file has just one contiguous domain and
# thus just one skeleton
for sk,skdist in me.skeletons:
    plot_3d(sk, color=skdist,figidx=1)



# some useful data, from the MorphExpt object
print "Radius of underlying droplet, obtained by fitting:", me.droplet[3]
print "Diameter of largest possible inscribed circle:", 2*np.max(me.sphere_distance())

# get the length of the longest piece of the backbone. see the skeleton_tools functions docs for more details.

gg = skt.skel_min_graph_real(me.skeletons[0],me.ptc_dia)
bb = skt.backbone(gg)
print "Length of backbone spine:", skt.skel_length(bb)


mlab.show()
