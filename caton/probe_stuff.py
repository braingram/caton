from __future__ import division, with_statement
import numpy as np, matplotlib.pyplot as plt, matplotlib.delaunay as sd
from string import strip
from utils_misc import index_safe, max2min, switch_ext
from utils_graphs import edges, nodes, complete_graph, add_edge, add_node
import re,string, matplotlib
from os.path import abspath, join, dirname, basename
from collections import namedtuple


PROBE_SITES = None
PROBE_DIM = None
PROBE_GRAPH = None
N_SITES = None
SORT_GROUPS = None
# Site class:
Site = namedtuple('Site', 'name dat x y')

def load_probe(ProbeFileName):
    """
    Loads site location and dat file mapping from probe file. 
    Each line of probe file has one of the following forms:
    sp1 0 (3.5 4.5)    # 2d location
    lp2 1 (4.5)    # 1d location
    tet1 3    # 0d location
    
    After this function is called, global variable PROBE_SITES contains a list of Site
    objects, where Site.name, Site.dat, Site.x, and Site.y hold name, dat-file channel,
    x-coor, y-coor (when relevant). The order of sites in the list is the order that the
    channels will appear after clustering."""
    #import wingdbstub

    global PROBE_SITES,PROBE_DIM,N_SITES,PROBE_GRAPH,SORT_GROUPS
        
    form0d = re.compile("^(\w+)\s(\d+)$")
    form1d = re.compile("^(\w+)\s(\d+)\s\(([-,\d,.]+)\)$")
    form2d = re.compile("^(\w+)\s(\d+)\s\(([-,\d,.]+)\s([-,\d,.]+)\)$")
    
    def read_form0d(line):
        match_obj = form0d.match(line)
        return Site(match_obj.group(1),int(match_obj.group(2)),None,None)
    def read_form1d(line):
        match_obj = form1d.match(line)
        return Site(match_obj.group(1),int(match_obj.group(2)),float(match_obj.group(3)),None)        
    def read_form2d(line):
        match_obj = form2d.match(line)
        return Site(match_obj.group(1),int(match_obj.group(2)),float(match_obj.group(3)),float(match_obj.group(4)))    
    
    def remove_comment(line):
        ind = line.find("#")
        return line if ind == -1 else line[:ind]
    
    with open(ProbeFileName,'r') as ProbeFile:
        print("Reading probe file: %s"%abspath(ProbeFileName))
        OrigLines = ProbeFile.readlines()

    Lines = map(remove_comment,OrigLines) #remove comments from lines
    Lines = map(strip,Lines) # strip leading and traiing whitespace
    Lines = filter(None,Lines) # get rid of empty lines
    
    edge_ind = index_safe("EDGES",Lines)
    group_ind = index_safe("GROUPS",Lines)
    SiteLines = Lines[:(edge_ind or group_ind)]
    
    if form0d.match(SiteLines[0]): 
        reader = read_form0d
        PROBE_DIM = 0
    elif form1d.match(SiteLines[0]): 
        reader = read_form1d
        PROBE_DIM = 1
    elif form2d.match(SiteLines[0]): 
        reader = read_form2d
        PROBE_DIM = 2
    else: 
        raise Exception("%s is not of valid form."%SiteLines[0])
        
    PROBE_SITES = map(reader,SiteLines)    
    N_SITES = len(PROBE_SITES)

    if edge_ind is None and group_ind is None:
        print("You didn't specify custom edges and groups in probe file, so I'm using defaults for %i-dimensional probe. If you want to customize edges and groups, start by copying and pasting the following text into your probe file."%PROBE_DIM)
        print("-----------------------")
        PROBE_GRAPH = make_graph()
        SORT_GROUPS = ch_subsets()
        
        print("".join(OrigLines))
        print("EDGES")
        print(string.join(get_edge_lines(),"\n"))
        print("GROUPS")
        print(string.join(get_group_lines(),"\n"))            
        print("-----------------------")        
    elif edge_ind is None or group_ind is None: raise Exception("Either specify both GROUPS and EDGES or neither.")
    else:
        print("Loading custom edges and groups from your probe file:")
        print("-----------------------")
        print("".join(OrigLines))        
        print("-----------------------")
        PROBE_GRAPH = edge_lines_to_graph(Lines[(edge_ind+1):group_ind])
        SORT_GROUPS = group_lines_to_groups(Lines[(group_ind+1):])

def edge_lines_to_graph(edge_lines):
    probe_graph = dict([(i_site,set([])) for i_site in xrange(N_SITES)])
    site_names = [site.name for site in PROBE_SITES]
    for line in edge_lines:
        src,targ = line.split()
        add_edge(probe_graph,site_names.index(src),site_names.index(targ))
    return probe_graph
def group_lines_to_groups(group_lines):
    site_names = [site.name for site in PROBE_SITES]    
    return [[site_names.index(name) for name in line.split()] for line in group_lines]
def get_edge_lines():  return ["%s %s"%(PROBE_SITES[src].name,PROBE_SITES[targ].name) for (src,targ) in edges(PROBE_GRAPH) if src<targ]
def get_group_lines(): return [string.join([PROBE_SITES[i_site].name for i_site in group]," ")
                               for group in SORT_GROUPS]
    
        
        
def ch_subsets():    
    if PROBE_DIM == 0 or N_SITES <= 4:
        return ch_subsets_0d()
    if PROBE_DIM == 1:
        ChLocs = np.array([site.x for site in PROBE_SITES])
        return ch_subsets_1d(ChLocs)
    if PROBE_DIM == 2:
        ChLocs = np.array([(site.x,site.y) for site in PROBE_SITES])
        return ch_subsets_2d(ChLocs)
    #if PROBE_DIM == "manual":
        #return SORT_GROUPS


def ch_subsets_0d():
    return [np.arange(N_SITES,dtype=np.int32)]
    
def ch_subsets_1d(ChLocsArr):
    SortInds = np.argsort(ChLocsArr)
    return [SortInds[i:(i+4)] for i in range(len(SortInds)-3)]

def ch_subsets_2d(ChLocsArr):
    T = sd.Triangulation(*ChLocsArr.transpose())
    triangles = T.triangle_nodes
    tri_adjacency = dict([(tri,nb_tris[nb_tris != -1]) for (tri,nb_tris) in enumerate(T.triangle_neighbors)])

    subsets = set()
    for src,targs in tri_adjacency.items():
        for targ in targs:
            if src < targ:
                subsets.add(tuple(set(triangles[src]).union(set(triangles[targ]))))
    return list(subsets)
        
        
        
#def load_probe_manual(Lines):
    #global PROBE_SITES,PROBE_DIM,N_SITES,SORT_GROUPS,PROBE_GRAPH
    #PROBE_SITES = []
    #PROBE_DIM = "manual"
    #PROBE_GRAPH = {}
    #SORT_GROUPS = []
    
    #for i_line,line in enumerate(Lines):
        #if line == "GROUPS": break
        
    #site_lines = Lines[:i_line]
    #group_lines = Lines[(i_line+1):]
    #print site_lines
    #print group_lines
    #names = [line.split()[0] for line in site_lines]
    
    #name2ind = lambda name: names.index(name)
    
    #for i_line,line in enumerate(site_lines):
        #words = line.split()
        #PROBE_SITES.append(Site(words[0],int(words[1]),None,None))
        #add_node(PROBE_GRAPH,name2ind(words[0]))
        #for neighbor in words[2:]:
            #add_edge(PROBE_GRAPH,name2ind(words[0]),name2ind(neighbor))
    #N_SITES = len(PROBE_SITES)
    
    #for line in group_lines:
        #SORT_GROUPS.append(map(lambda name: name2ind(name),line.split()))
    
def make_graph():
    if PROBE_DIM == 0:
        return complete_graph(N_SITES)
    elif PROBE_DIM == 1:
        x = [site.x for site in PROBE_SITES]
        return path_graph(x)
    elif PROBE_DIM == 2:
        xy = [(site.x,site.y) for site in PROBE_SITES]
        return triangle_graph(xy)


def plot_probe(ProbeFileName,output_dir=None):

    load_probe(ProbeFileName)
    
    if PROBE_DIM == 0:
        x = [np.cos(2*np.pi*i_site/N_SITES) for i_site in xrange(N_SITES)]
        y = [np.sin(2*np.pi*i_site/N_SITES) for i_site in xrange(N_SITES)] 
    elif PROBE_DIM == 1:
        x = [site.x for site in PROBE_SITES]
        y = [0 for site in PROBE_SITES]
    elif PROBE_DIM == 2:
        x = [site.x for site in PROBE_SITES]
        y = [site.y for site in PROBE_SITES]

    #if max2min(x) > max2min(y):
        #figsize = (6,1+max2min(y)/max2min(x)*6)
    #else:
        #figsize = (1+max2min(x)/max2min(y)*6,6)        
               
    plt.plot(x,y,'go')
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    
    for ind_src,ind_targ in edges(PROBE_GRAPH):
        ax.add_line(matplotlib.lines.Line2D([x[ind_src],x[ind_targ]],[y[ind_src],y[ind_targ]],lw=1))
    
    for ind_node in nodes(PROBE_GRAPH):
        site = PROBE_SITES[ind_node]
        ax.text(x[ind_node],y[ind_node],site.name,color='r')
            
    
    img_filename = join(output_dir or dirname(ProbeFileName),switch_ext(basename(ProbeFileName),"png"))
    print("Saving figure as %s"%abspath(img_filename))
    plt.savefig(img_filename)

    
def triangle_graph(Locs):
    """
    Returns a graph giving the Delaunay triangulation of a set of two-dimensional points    

    Parameters
    -------
    Locs : ndarray, or anything that gets cast into ndarray upon np.array(Locs)
        Locs.shape = (n,2), where n is the number of points

    Returns
    ------
    out : networkx Graph
    """
    Locs = np.array(Locs)
    if len(Locs)==1:
        return complete_graph(1)
    else:
        Triangulation = sd.Triangulation(Locs[:,0],Locs[:,1])
        return Triangulation.node_graph()

    
def path_graph(Locs):
    """
    Returns a line graph for a set of one-dimensional locations.    

    Parameters
    -------
    Locs : ndarray, or anything that gets cast into ndarray upon np.array(Locs)

    Returns
    ------
    out : graph, i.e. dictionary key -> set(keys)
    """
    Locs = np.array(Locs).flatten()
    if len(Locs)==1:
        return complete_graph(1)    
    else:
        G = {}
        SortInds = np.argsort(Locs)
        for src,targ in zip(SortInds[1:],SortInds[:-1]):
            add_edge(G,src,targ)
        return G