from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

INT32 = np.int32
ctypedef np.int32_t INT32_t
FLOAT32 = np.float32
ctypedef np.float32_t FLOAT32_t
INT8 = np.int8
ctypedef np.int8_t INT8_t
FLOAT64 = np.float64
ctypedef np.float64_t FLOAT64_t

cdef float HUGE = 1e30
#cdef extern from "math.h":
#    double log(double)
log = np.log
cdef float pi = 3.141592

def transitions(np.ndarray[INT32_t,ndim=1] X_n,np.ndarray[INT32_t,ndim=1] Y_n,int M):

    cdef np.ndarray[INT32_t,ndim=2] transitions = np.zeros((M,M),INT32)
    cdef int m1,m2,N
    N = X_n.size
    for n in range(N):
        transitions[X_n[n],Y_n[n]] += 1
    return transitions
    
def sqnorm(np.ndarray[FLOAT32_t,ndim=2] Arr_mn):
    """Map each coumn onto squared norm"""
    cdef int n,m,N,M
    M = Arr_mn.shape[0]
    N = Arr_mn.shape[1]
    cdef np.ndarray[FLOAT32_t,ndim=1] Out = np.zeros(N,FLOAT32)  

    for m in range(M):
        for n in range(N):
            Out[n] += Arr_mn[m,n]**2

    return Out

def arr_for_c(graph_dic):
    n_nodes = len(graph_dic)
    assert graph_dic.keys() == range(n_nodes)
    edges_arr = np.zeros((n_nodes,n_nodes),dtype=np.int32)
    n_edges_arr = np.zeros(n_nodes,dtype = np.int32)
    for (source,targs) in graph_dic.items():
        n_edges_arr[source] = len(targs)+1 # +1 including itself
        edges_arr[source,0:n_edges_arr[source]] = [source]+list(targs)
    return edges_arr,n_edges_arr

def argmax_two_per_row(np.ndarray[FLOAT32_t,ndim=2] Arr_nm):
    cdef int M,N,m,n
    N,M = Arr_nm.shape[0],Arr_nm.shape[1]
    cdef np.ndarray[INT32_t,ndim=1] First = np.empty(N,INT32)
    cdef np.ndarray[INT32_t,ndim=1] Second = np.empty(N,INT32)    
    cdef int arg_first, arg_second, arg_new
    cdef float first, second, new    
    
    for n in range(N):
        first,second = -np.inf,-np.inf
        argfirst,argsecond = 0,0
        for m in range(M):
            new,argnew = Arr_nm[n,m], m
            if new > second:
                if new > first:
                    first,second = new,first
                    argfirst,argsecond = argnew,argfirst
                else:
                    second = new
                    argsecond = argnew
        First[n],Second[n] = argfirst,argsecond        
    return First,Second

    
def class_means(np.ndarray[FLOAT32_t,ndim=2] X_nf, np.ndarray[INT32_t,ndim=1] Class_n, int M):
    cdef int F,N,f,m,n
    N,F = X_nf.shape[0],X_nf.shape[1]
    cdef np.ndarray[FLOAT64_t,ndim=2] Mu_mf = np.zeros((M,F),FLOAT64)
    cdef np.ndarray[INT32_t,ndim=1] ClassCounts_m = bincount(Class_n,M)
    
    for n in range(N):
        for f in range(F):
            Mu_mf[Class_n[n],f] += X_nf[n,f]
    
    for m in range(M):
        if ClassCounts_m[m] > 0:
            Mu_mf[m] /= ClassCounts_m[m]
    
    return Mu_mf.astype(FLOAT32)

    

def class_covs(np.ndarray[FLOAT32_t,ndim=2] X_nf, np.ndarray[FLOAT32_t,ndim=2] Mean_mf, np.ndarray[INT32_t,ndim=1] Class_n, int M):
    cdef int F,N,f1,f2,m,n
    N,F = X_nf.shape[0],X_nf.shape[1]
    cdef np.ndarray[FLOAT64_t,ndim=3] Cov_mff = np.tile(.01*np.eye(F,dtype=FLOAT64),(M,1,1))
    cdef np.ndarray[FLOAT64_t,ndim=1] Vec2Mean
    ClassCounts_m = bincount(Class_n,M)    
    
    for n in range(N):
        Vec2Mean = (X_nf[n] - Mean_mf[Class_n[n]]).astype(FLOAT64)
        for f1 in range(F):
            for f2 in range(F):
                Cov_mff[Class_n[n],f1,f2] += Vec2Mean[f1]*Vec2Mean[f2]
    
    for m in range(M):
        if ClassCounts_m[m] > 0:
            Cov_mff[m] /= ClassCounts_m[m]
    
    return Cov_mff.astype(FLOAT32)

    
def class_wts(np.ndarray[INT32_t,ndim=1] Class_n,M):
    return bincount(Class_n,M).astype(np.float32) / len(Class_n)
    
    
def class_sums(np.ndarray[FLOAT32_t,ndim=1] Y_n, np.ndarray[INT32_t,ndim=1] Class_n, M):
    cdef int N = Y_n.size
    cdef np.ndarray[FLOAT32_t,ndim=1] Sum_m = np.zeros(M,FLOAT32)    
    for n in range(N):
        Sum_m[Class_n[n]] +=  Y_n[n] 
    return Sum_m
    

def bincount(np.ndarray[INT32_t,ndim=1] Class_n, M):
    cdef np.ndarray[INT32_t,ndim=1] Count_m = np.zeros(M,INT32)
    cdef int n
    for n in range(len(Class_n)):
        Count_m[Class_n[n]] += 1
    return Count_m
    
def subset_inds(Class_n,int M=0):
    M = M or Class_n.max()+1
    IndLists = [[] for i in range(M)]
    for i,cl in enumerate(Class_n):
        IndLists[cl].append(i)
    return IndLists

    
def connected_components(np.ndarray[INT8_t,ndim=2] st_arr, ch_graph, int s_back):
    """ blah blah blah"""
    cdef int n_s = st_arr.shape[0]
    cdef int n_ch = st_arr.shape[1]
    cdef int n_s_buff = s_back+1
    cdef np.ndarray[INT32_t, ndim=2] label_buffer = np.zeros((n_s_buff,n_ch),dtype=INT32)  
    cdef int i_s, i2_s, i_ch, j_s, j_ch, j_ch_ind, j_sstart,
    cdef INT32_t adjlabel,c_label
    cdef dict comp_inds = {}
    
    edges_arr_temp,n_edges_arr_temp = arr_for_c(ch_graph)
    cdef np.ndarray[INT32_t, ndim=2] edges_arr = edges_arr_temp.copy()
    cdef np.ndarray[INT32_t, ndim=1] n_edges_arr = n_edges_arr_temp.copy()    
    
    c_label = 1    
    
    for i_s in range(n_s):

        i2_s = i_s % n_s_buff        
        for i_ch in range(n_ch):
            label_buffer[i2_s,i_ch] = 0
        
        for i_ch in range(n_ch):
            if st_arr[i_s,i_ch]:

                for j_s in range(0,n_s_buff):
                    for j_ch_ind in range(n_edges_arr[i_ch]):
                        j_ch = edges_arr[i_ch,j_ch_ind]
                        
                        if label_buffer[j_s,j_ch] != 0: # adjacent element is nonzero
                            adjlabel = label_buffer[j_s,j_ch]
                            if label_buffer[i2_s,i_ch] == 0: # current element is still zero
                                label_buffer[i2_s,i_ch] = adjlabel
                                comp_inds[adjlabel].append((i_s,i_ch))
                            elif label_buffer[i2_s,i_ch] != adjlabel:
                                samps_chans = np.array(comp_inds[adjlabel])
                                samps_chans = samps_chans[i_s - samps_chans[:,0] < n_s_buff]
                                label_buffer[samps_chans[:,0]%n_s_buff,samps_chans[:,1]] = label_buffer[i2_s,i_ch]
                                comp_inds[label_buffer[i2_s,i_ch]].extend(comp_inds.pop(adjlabel))
                    # nothing adjacent
                if label_buffer[i2_s,i_ch] == 0:
                    label_buffer[i2_s,i_ch] = c_label
                    comp_inds[c_label] = [(i_s,i_ch)]                    
                    c_label += 1

    return comp_inds.values()      