#!/usr/bin/env python
from __future__ import division
import numpy as np
from optparse import OptionParser
from scipy.signal import lfilter,butter,lfilter_zi
from caton.progressbar import ProgressBar
from caton.utils_misc import switch_ext
from caton.core import get_dat_pars
from os.path import abspath

# Defaults
DT_IN = "i2"
DT_OUT = "i2"
RATIO = 16
BUTTER_ORD = 3

def filt(fname_in,fname_out,two_way,f_hi=None,f_lo=None):
    print "Input file: %s"%fname_in
    print "Output file: %s"%fname_out
    
    n_ch,sample_rate = get_dat_pars(fname_in)        
    n_ch,_ = get_dat_pars(fname_in)    
    in_memmap = np.memmap(fname_in,mode="r",dtype=DT_IN).reshape(-1,n_ch)
    n_s = in_memmap.shape[0] 
    out_memmap = np.memmap(fname_out,mode="w+",shape=(n_s,n_ch),dtype=DT_OUT)
    
    if f_lo is None: b,a = butter(BUTTER_ORD,f_hi/(sample_rate/2),"high")
    elif f_hi is None: b,a = butter(BUTTER_ORD,f_hi/(sample_rate/2),"low")
    else: b,a = butter(BUTTER_ORD,(f_lo/(sample_rate/2),f_hi/(sample_rate/2)),"pass")

    if two_way: filter_2way(b,a,in_memmap, out_memmap)
    else: filter_1way(b,a,in_memmap,out_memmap)

def filter_1way(b,a,in_memmap,out_memmap):
    
    X_nc = in_memmap
    n_s, n_ch = X_nc.shape
    Y_mc = out_memmap
    assert Y_mc.shape == (n_s,n_ch)
    
    n_out_chunk = 10000
    n_in_chunk = n_out_chunk

    t_warmup = max(a.size,b.size)*3
    
    pbar = ProgressBar(maxval = n_s).start()
    z = np.zeros((max(a.size,b.size)-1,n_ch))
    for i_ch in xrange(n_ch):
        _,z[:,i_ch] = lfilter(b,a,X_nc[t_warmup-1::-1,i_ch],zi=lfilter_zi(b,a))
    for i_in,i_out in zip(xrange(0,n_s,n_in_chunk),xrange(0,n_s,n_in_chunk)):
        for i_ch in xrange(n_ch):
            Y_mc[i_out:i_out+n_out_chunk,i_ch],z[:,i_ch] = lfilter(b,a,X_nc[i_in:i_in+n_in_chunk,i_ch],zi=z[:,i_ch])
        pbar.update(i_in)
    pbar.finish()

def filter_2way(b,a,in_memmap,out_memmap):
    print "Forward filtering"
    filter_1way(b,a,in_memmap,out_memmap)
    print "Backward filtering"
    filter_1way(b,a,out_memmap[::-1],out_memmap[::-1])

    

if __name__ == '__main__':
    usage = """
    Filter your data between f_hi and f_lo, measured in Hz (the frequencies you give will be compared with sample rate in xml file.)
    Give either f_hi, f_lo, or both using --hi and --lo.
    By default this script does a two-way filter (like filtfilt).
    If output file is not specified, it will have extension fil (foo.dat -> foo.fil)
    
    %prog your_dat_file.dat [output filename] [--hi=f_hi] [--lo=f_lo]
    %prog -h displays help"""
    parser = OptionParser(usage)
    parser.add_option("-1",help="One-way filter.",action="store_false",default=True,dest="two_way")
    parser.add_option("--hi",help="Top of passband (Hz)",type=float,action="store",dest="hi")
    parser.add_option("--lo",help="Bottom of passband (Hz)",type=float,action="store",dest="lo")
    parser.add_option("--in-type",help="Datatype of input file. One of i2,i4,i8,f4,f8. Default i2 (16-bit integer).",action="store",dest="dt_in",default=DT_IN)
    parser.add_option("--out-type",help="Datatype of output file. One of i2,i4,i8,f4,f8. Default i2 (16-bit integer).",action="store",dest="dt_out",default=DT_OUT)
    (opts,args) = parser.parse_args()
    if opts.hi is None and opts.lo is None:
        parser.error("must specify high pass freq or low pass freq with --hi or --lo")
    OUTNAME = args[1] if len(args) > 1 else switch_ext(abspath(args[0]),"fil")
    DT_IN = opts.dt_in
    DT_OUT = opts.dt_out
    
    filt(args[0],OUTNAME,two_way=opts.two_way,f_hi=opts.hi,f_lo=opts.lo)