#!/usr/bin/env python

from caton.core import get_dat_pars
import caton.probe_stuff as probe_stuff
import matplotlib.pyplot as plt, numpy as np
import os, caton.myopts
from scipy.signal import lfilter,butter
from optparse import OptionParser
from caton.utils_misc import find_file_with_ext,basename_noext
from caton.utils_graphs import edges

TINY = 1e-4

def zerodiag(arr):
    out = arr.copy()
    out[xrange(arr.shape[0]),xrange(arr.shape[0])] = 0
    return out

def check_crosstalk(DatFileName,ProbeFileName=None,output_dir = None):
    "Makes a colormap that shows conditional correlations between channels. Output it as a png file."
    N_SAMPLES = 80000
    if not os.path.exists(DatFileName):
        raise Exception("Dat file %s does not exist"%DatFileName)
    DatFileName = os.path.abspath(DatFileName)       
    DatDir = os.path.abspath(os.path.dirname(DatFileName))
    n_ch_dat,sample_rate = get_dat_pars(DatFileName)

    ProbeFileName = ProbeFileName or find_file_with_ext(DatDir,'probe',ex_if_not_found=False)
    if ProbeFileName is not None:
        probe_stuff.load_probe(ProbeFileName)        
        Channels_dat = [site.dat for site in probe_stuff.PROBE_SITES]
        ChNames = [site.name for site in probe_stuff.PROBE_SITES]
        ChannelGraph = probe_stuff.PROBE_GRAPH            
    else:
        Channels_dat = range(n_ch_dat)
        ChNames = [str(i) for i in range(n_ch_dat)]

    
        
    b,a = butter(3,(600./(sample_rate/2),.95),'pass')
    X_nc = np.fromfile(DatFileName,dtype=np.int16,count=N_SAMPLES*n_ch_dat).reshape(-1,n_ch_dat)[:,Channels_dat]
    X_nc = lfilter(b,a,X_nc,axis=0)
    Cov_cc = np.cov(X_nc.T)
    try: Con_cc = np.linalg.inv(Cov_cc)
    except np.linalg.LinAlgError:
        vals,vec_mat = np.linalg.eigh(Cov_cc)
        print "Singular covariance matrix. Small eigvals/eigvecs:"
        np.set_printoptions(precision=2,suppress=True)
        for _val,vec in zip(filter(lambda x: x < TINY, vals),vec_mat.T):
            print vec
        exit()
    
    n_ch = len(Channels_dat)                
    plt.pcolor(np.ma.array(np.abs(Con_cc),mask=np.eye(Con_cc.shape[0]))
                           ,cmap=plt.cm.cool)
    #plt.colorbar()
    
    
    if ProbeFileName is not None and len(edges(ChannelGraph))>0:
        ptsx1,ptsy1 = np.array(edges(ChannelGraph)).T
        plt.plot(ptsx1+.5,ptsy1+.5,'kx',ms=10,mew=3)
        
        
    ax = plt.gca()
    ax.set_xticks(np.arange(n_ch)+.5)
    ax.set_yticks(np.arange(n_ch)+.5)
    ax.set_xticklabels(ChNames,rotation='vertical')
    ax.set_yticklabels(ChNames)
    
    basename = basename_noext(ProbeFileName) if ProbeFileName else basename_noext(DatFileName)
    outpath = os.path.join(output_dir or ProbeFileName and os.path.dirname(ProbeFileName) or DatDir,"%s_crosstalk.png"%basename)
    print("saving figure in %s"%outpath)
    plt.savefig(outpath)   
    
if __name__ == '__main__':
    usage = """%prog your_dat_file.dat [options]
    %prog -h displays help"""
    parser = OptionParser(usage)
    parser.add_options([caton.myopts.probe,caton.myopts.output])
    (opts,args) = parser.parse_args()
    check_crosstalk(args[0],opts.probe,opts.output)
