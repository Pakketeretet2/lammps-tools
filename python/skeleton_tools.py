import sys, numpy
import scipy
from scipy import optimize
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

from networkx import nx
from copy import copy,deepcopy
from numpy import newaxis


hull_path = "./hull"


def droplet_fitting(coords):
    '''
    find the droplet substrate from the coordinates of particles:
    coords: coordinates of partilces on the droplet interface, in [x, y, z]
    return the droplet size (radius) to the center of particles.
    '''
    print >> sys.stderr, "Someone called droplet_fitting..."
    # fitting function, calculate the distance of a point (pos[x, y, z])
    # from the sphere (center at s[:3], with radius s[3]).
    s0=[0, 0, 0, 10] #initial guess of the sphere center and radius
    if coords.shape[0]<4:
        return s0
    fitfunct=lambda s, x, y, z:(x-s[0])**2+(y-s[1])**2+(z-s[2])**2-s[3]**2
    s, success=optimize.leastsq(fitfunct, s0, args=(coords[:,0], coords[:,1], coords[:,2]))
    s[3]=abs(s[3])
    return s

    

def fit_radii(coords,droplet):
    coords = coords - droplet[:3]
    return droplet, numpy.sqrt(numpy.sum(coords**2,axis=1))

    
def binarify(arr):
    arr[numpy.where(arr > 0)] = 1
    return arr
    
def toPoints(arr):
    return numpy.transpose(numpy.where(arr));
    
# visualization helpers
def get_central_vectors(pts):
    drp = droplet_fitting(pts)
    pts = pts-drp[:3]
    u,v,w = pts.T
    return u,v,w
    
    
def detect_local_minima(arr):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    background = (arr==0)
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    detected_minima = local_min - eroded_background
    return numpy.where(detected_minima)
    
def detect_local_maxima(arr):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    local_max = (filters.maximum_filter(arr, footprint=neighborhood)==arr)
    background = (arr==0)
    # 
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    detected_maxima = local_max - eroded_background
    return detected_maxima


def eliminate_rays_dist(arr,distmap,delta,gridpoints=[0],drp=None):
    """shoot rays from center of best fit sphere to points and eliminate
    those that overlap in some region of theta, phi"""
    binarr = deepcopy(arr)
    pts = numpy.transpose(numpy.where(binarr))
    
    if len(pts) < 6: return binarr, pts, distmap[numpy.where(binarr)] # if too small, fitting sphere is bad
    
    if drp is None:
        if numpy.rank(gridpoints) == 3:
            drp = droplet_fitting(numpy.transpose(numpy.where(gridpoints)))
        else:
            drp = droplet_fitting(pts)
    rad = drp[3]
    pts0 = pts-drp[:3]
    x,y,z = numpy.transpose(pts0)
    r = numpy.sqrt(x**2+y**2+z**2)
    theta = numpy.arccos(z/r)
    phi = numpy.arctan2(y,x)
    
    dropdist = numpy.abs(r-rad)
    order = dropdist.argsort()
    
    # remove points that are too close to others in spherical coords
    for idx in order:
        pt = pts[idx]
        if binarr[tuple(pt)]:
            distances = []
            th = theta[idx]
            ph = phi[idx]
            delth = delta/rad/2.
            delph = delta/(rad*numpy.sin(th))/2.
            #~ delth = delta/rad/2.
            boxed = (th - delth < theta) & (theta < th +delth) & (ph-delph < phi) & (phi < ph+delph)
            tempdist = []
            for i in numpy.where(boxed)[0]:
                tempdist.append(distmap[tuple(pts[i])])
                if i != idx: binarr[tuple(pts[i])] = 0
            distmap[tuple(pt)] = numpy.max(tempdist) # assign maximum distance to the kept pt
    
    # generate projected points and return as well
    ptsprime = numpy.transpose(numpy.where(binarr))
    ptsprime0 = ptsprime-drp[:3]
    x,y,z = numpy.transpose(ptsprime0)
    r = numpy.sqrt(x**2+y**2+z**2)
    theta = numpy.arccos(z/r)
    phi = numpy.arctan2(y,x)
    
    xs = rad * numpy.cos(phi)*numpy.sin(theta)
    ys = rad * numpy.sin(phi)*numpy.sin(theta)
    zs = rad * numpy.cos(theta)
    
    return binarr, numpy.transpose(numpy.array([xs,ys,zs]))+drp[:3], distmap[numpy.where(binarr)]
    #~ return binarr, ptsprime, distmap[numpy.where(binarr)]


def get_skeletons_binary(eucb, grid_domains, resolution_scale,ptc_dia):
    ''' 
    return a tuple of binary arrays, each identifying the regions of
    large negative laplacian in eucb, with some size/connectivity requirements
    '''
    lap = filters.laplace(eucb)
    lap[eucb < (resolution_scale*ptc_dia)] = 0 # exclude points near the bdry
    lap[lap > -0.8] = 0 # medial axis is the large negative laplacian areas (value -1 because the length scale of variation is one lattice site) (added a tolerance so slightly more than -1 ok)
    
    
    skell = numpy.zeros_like(lap).astype(int)
    skell[(lap < 0) & (grid_domains > 0)] = 1
    
    skell3 = morphology.binary_dilation(skell)
    
    res1, nel = scipy.ndimage.measurements.label(skell3)
    sizes = [numpy.sum(res1 == i) for i in numpy.arange(1,nel+1)]
    
    skeleton_closet = [] # collect binary arrays of distinct skeletons
    for i in numpy.arange(1,nel+1):
        size = sizes[i-1]
        if (size > (6*ptc_dia*resolution_scale*3)):
            newskel = numpy.zeros_like(skell3)
            newskel[numpy.where(res1 == i)] = 1
            skeleton_closet.append(morphology.binary_erosion(newskel))
            
    if len(skeleton_closet) > 0: return tuple(skeleton_closet)
    else:
        #~ idx = numpy.argmax(sizes) # if no large skeleton found, return the largest one, can't do better
        skellprime = numpy.zeros_like(lap).astype(int)
        skellprime[numpy.unravel_index(numpy.argmin(lap),lap.shape)] = 1
        return tuple([skellprime])
    

def get_skeletons(eucb, grid_domains,resolution_scale,ptc_dia,drp=None):
    """
    return a tuple of skeletons. Each skeleton itself is a tuple, the 
    first element being a set of points, and the second element the
    corresponding distance-transform value for that point on the skeleton.
    """
    delta = resolution_scale*ptc_dia/2.1
    skels = get_skeletons_binary(eucb, grid_domains,resolution_scale,ptc_dia)
    output = []
    for skel in skels:
        res, respts, rds = eliminate_rays_dist(skel,eucb,delta=delta, gridpoints=grid_domains, drp=drp)
        output.append((respts,rds))
    return tuple(output)
    
    
def skel_min_graph(skeleton,resolution_scale,ptc_dia,tol=0.5):
    """
    given a skeleton, return a minimum spanning tree that connects the
    points in the skeleton using an adjacency criterion. the weights of
    the edges in the graph correspond to real distances between the points
    """
    delta = resolution_scale*ptc_dia/1.8
    coords, rads = skeleton
    N = coords.shape[0]
    
    adjm = numpy.zeros([N, N], dtype=int)
    dist_mat = numpy.array((coords[0:,newaxis] - coords))
    dist_mat = numpy.sqrt(numpy.sum(dist_mat**2,axis=2))
    adjm = (dist_mat < (1+tol)*delta)
    adjm = adjm - numpy.eye(N)
    adjm = adjm.astype('bool')
    
    g = nx.Graph(dist_mat*adjm)
    return nx.minimum_spanning_tree(g)
    
def skel_min_graph_real(skeleton,ptc_dia,tol=0.5):
    """
    same as skel_min_graph, but using real coordinates and distances
    as opposed to lattice
    """
    delta = ptc_dia/1.8
    coords, rads = skeleton
    N = coords.shape[0]
    
    adjm = numpy.zeros([N, N], dtype=int)
    dist_mat = numpy.array((coords[0:,newaxis] - coords))
    dist_mat = numpy.sqrt(numpy.sum(dist_mat**2,axis=2))
    adjm = (dist_mat < (1+tol)*delta)
    adjm = adjm - numpy.eye(N)
    adjm = adjm.astype('bool')
    
    g = nx.Graph(dist_mat*adjm)
    for i in range(len(rads)):
        g.node[i]['rad'] = rads[i]
    if not nx.is_connected(g): # return largest subgraph if disconnected
        # g = nx.connected_component_subgraphs(g)[0]
        g = sorted(nx.connected_component_subgraphs(g), key = len, reverse=True)[0]
    return nx.minimum_spanning_tree(g)


def condense_skeletons(skeletons):
    ''' condense all skeletons to get one skeleton of the domain'''
    sks = [sk[0] for sk in skeletons]
    skdists = [sk[1] for sk in skeletons]
    return (numpy.concatenate(sks),numpy.concatenate(skdists))

def find_endpts(g):
    """
    return node indices for nodes of degree one
    """
    return [node for node in g.nodes_iter() if (g.degree(node) == 1)]

def find_branches(g):
    """
    return the endpoints, corresponding junction nodes, and length to
    the junction nodes for the graph g
    """
    juncts = []
    branchlengths = []
    endpts = find_endpts(g)
    for endpt in endpts:
        junct = endpt
        prev = endpt
        branchlength = 0
        for j in nx.dfs_preorder_nodes(g,endpt):
            if j != endpt:
                branchlength += g[j][prev]['weight']
                if g.degree(j) > 2:
                    junct = j
                    break
            prev = j
        juncts.append(junct)
        branchlengths.append(branchlength)
    return endpts,juncts, branchlengths
    
def prune_smallest_branch(g):
    """
    remove the nodes and edges corresponding to the smallest branch
    by length (summed weights of the edges)
    """
    endpts, juncts, branchlengths = find_branches(g)
    if len(endpts) == 2: return 0
    if len(branchlengths) == 0: return 0    
    smidx = numpy.argmin(branchlengths)
    endpt = endpts[smidx]
    junct = juncts[smidx]
    branch = nx.shortest_path(g,endpt,junct)
    g.remove_nodes_from(branch[:-1])
    return 1
    
def prune_all_branches(g):
    """
    prune branches until only the longest branch --  the "backbone" --- remains
    """
    while prune_smallest_branch(g): continue
    
def backbone(g):
    g2 = deepcopy(g)
    prune_all_branches(g2)
    return g2
        
        
# the following functions should only be run on 
def skel_length(g):
    """
    add up all bond lengths within graph. if you want this to mean a meaningful length,
    you must have already pruned all the edges in g
    """
    if len(find_endpts(g)) > 2:
        print "Warning: running skel_length() on a skeleton with branches"
    return numpy.sum([g[i][j]['weight'] for i,j in g.edges()])
    
def get_tangents(g,pts):
    """
    return a set of tangent vectors connecting successive points
    """
    ep = find_endpts(g)
    if len(ep) > 2:
        print "Warning: running get_angles() on a skeleton with branches"
    
    bb = nx.dfs_preorder_nodes(g, ep[1])
    p1 = bb.next()
    tangents = []
    for p2 in bb:
        tangents.append(pts[p2]-pts[p1])
        p1 = p2
    return tangents
    
def get_unit_tangents(g,pts):
    """
    return a set of unit tangent vectors connecting successive points
    """
    tg = get_tangents(g,pts)
    return [t/numpy.sqrt(numpy.dot(t,t)) for t in tg]
    
def get_lengths(g,pts):
    """
    return the length of each tangent
    """
    return [numpy.sqrt(numpy.dot(t,t)) for t in get_tangents(g,pts)]

def get_angles(g,pts):
    """
    return all the angles between pairs of successive segments of the graph.
    """
    def angle(p1,p2,p3):
        """return the angle between v32 and v21, tangent-tangent"""
        v23 = pts[p3]-pts[p2]
        v12 = pts[p2]-pts[p1]
        return numpy.arccos(numpy.dot(v12,v23)/numpy.sqrt(numpy.dot(v12,v12)*numpy.dot(v23,v23)))
    
    ep = find_endpts(g)
    if len(ep) > 2:
        print "Warning: running get_angles() on a skeleton with branches"
    
    bb = nx.dfs_preorder_nodes(g, ep[0])
    p1 = bb.next()
    p2 = bb.next()
    angles = []
    for p3 in bb:
        angles.append(angle(p1,p2,p3))
        p1 = p2
        p2 = p3
    return numpy.array(angles)
        
def ttcorr(g,pts):
    """
    return the tangent-tangent correlation function
    """
    tgs = get_tangents(g,pts)
    lens = [numpy.sqrt(numpy.dot(t,t)) for t in tgs]
    unittgs = numpy.array([t/numpy.sqrt(numpy.dot(t,t)) for t in tgs])
    
    dpmx = numpy.sum(unittgs[:,newaxis]*unittgs,axis=2)
    ttshift = [numpy.roll(dpmx[i],-i) for i in range(len(dpmx))]
    # pad with zeros?
    for i in range(1,len(ttshift)):
        ttshift[i][-i:] = 0
    corrfn = numpy.sum(ttshift,axis=0)
    return corrfn/len(corrfn)
    #~ return corrfn/range(len(corrfn),0,-1) # compensate for padding?
        
    
# various aspect ratios
def ar_backbone(g):
    """
    return aspect ratio defined by the length of the backbone (result of
    pruning all branches of g) divided by the diameter of the largest 
    incircle along the branch (twice the largest value of the distance
    transforms on the backbone)
    """
    dt = nx.get_node_attributes(g,'rad').values()
    diameter = 2*numpy.max(dt)
    return skel_length(backbone(g))/diameter
    
def dia_backbone_average(g):
    """
    return average diameter of the incircles along the longest
    branch
    """
    radii = numpy.array(nx.get_node_attributes(g,'rad').values())
    ave_dia = numpy.average(2*radii)
    return ave_dia    
    
def ar_backbone_average(g):
    """
    return aspect ratio defined by the length of the backbone
    divided by the average diameter of the incircles along the longest
    branch
    """
    return skel_length(backbone(g))/dia_backbone_average(g)

def dia_backbone_weighted(g):
    """
    return *weighted average* diameter of the incircles along the longest
    branch weighted by the area of each incircle
    """
    radii = numpy.array(nx.get_node_attributes(g,'rad').values())
    ave_dia = numpy.average(2*radii, weights = radii**2)
    return ave_dia

def ar_backbone_weighted(g):
    """
    return aspect ratio defined by the length of the backbone
    divided by the *weighted average* diameter of the incircles along the longest
    branch weighted by the area of each incircle
    """
    return skel_length(backbone(g))/dia_backbone_weighted(g)
    

    
