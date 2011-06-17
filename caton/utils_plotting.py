from __future__ import division
import numpy as np
from utils_misc import *
import matplotlib.pyplot as plt

def plot2d(X_sc,Ax):
    p2p = max2min(X_sc.flatten())
    n_s,n_ch = X_sc.shape
    for i_ch in xrange(n_ch):        
        Ax.plot(X_sc[:,i_ch]-i_ch*p2p/2)
        
def ccg_hist(Tm_n,sample_rate):
    t_diff = .050 # 50 miliseconds
    s_diff = int(sample_rate*t_diff)
    s_max = Tm_n.max()
    Binary_n = np.zeros(s_max+2*s_diff-1,dtype=np.int32)
    pos2signed = np.arange(-s_diff,s_diff,dtype=np.int32)
    Binary_n[Tm_n+s_diff] = 1
    hist = np.correlate(Binary_n,Binary_n[s_diff:-s_diff+1],mode='valid')
    hist[pos2signed.tolist().index(0)] = 0
    tseries = pos2signed / sample_rate
    return hist,tseries
    
def plot_ccg(Tm_n,Ax,sample_rate = 20000):
    hist,tseries = ccg_hist(Tm_n,sample_rate)
    Ax.plot(tseries,hist)
    
def test_plot2d():
    t = np.arange(0,10,.1)
    x = np.sin(np.vstack((t,t,t))).T
    Ax = plt.figure().add_subplot(1,1,1)
    plot2d(np.sin(x),Ax)
    plt.show()