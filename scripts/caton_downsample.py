#!/usr/bin/env python
from __future__ import division
import numpy as np
from optparse import OptionParser
from scipy.signal import lfilter,butter,lfilter_zi
from caton.progressbar import ProgressBar
from caton.utils_misc import switch_ext
from caton.core import get_dat_pars

# Defaults
DT_IN = "i2"
DT_OUT = "i2"
RATIO = 16
BUTTER_ORD = 3

def downsample_1way(in_memmap,out_memmap,ratio=16):
    b,a = butter(BUTTER_ORD,1/(2*ratio),'low')
    
    X_nc = in_memmap
    n_s, n_ch = X_nc.shape
    n_s2 = n_s//ratio
    Y_mc = out_memmap
    assert Y_mc.shape == (n_s2,n_ch)
    
    n_out_chunk = 10000
    n_in_chunk = n_out_chunk*ratio

    t_warmup = max(a.size,b.size)*3
    
    pbar = ProgressBar(maxval = n_s).start()
    z = np.zeros((BUTTER_ORD,n_ch))
    for i_ch in xrange(n_ch):
        _,z[:,i_ch] = lfilter(b,a,X_nc[t_warmup-1::-1,i_ch],zi=lfilter_zi(b,a))
    for (i_in,i_out) in zip(xrange(0,n_s,n_in_chunk),xrange(0,n_s2,n_out_chunk)):
        for i_ch in xrange(n_ch):
            Low_nc,z[:,i_ch] = lfilter(b,a,X_nc[i_in:i_in+n_in_chunk,i_ch],zi=z[:,i_ch])
            Y_mc[i_out:i_out+n_out_chunk,i_ch] = Low_nc[::ratio]
        pbar.update(i_in)
    pbar.finish()

def downsample_2way(in_memmap,out_memmap,ratio=16):
    print "Forward filtering"
    downsample_1way(in_memmap,out_memmap,ratio=ratio)
    print "Backward filtering"
    downsample_1way(out_memmap[::-1],out_memmap[::-1],ratio=1)
    
def downsample(fname_in,fname_out,two_way, ratio = 16):
    print "Input file: %s"%fname_in
    print "Output file: %s"%fname_out
    n_ch,_ = get_dat_pars(fname_in)    
    in_memmap = np.memmap(fname_in,mode="r",dtype=DT_IN).reshape(-1,n_ch)
    n_s = in_memmap.shape[0] 
    out_memmap = np.memmap(fname_out,mode="w+",shape=(n_s//ratio,n_ch),dtype=DT_OUT)
    if two_way: downsample_2way(in_memmap, out_memmap,ratio=16)
    else: downsample_1way(in_memmap,out_memmap,ratio=16)
    

def test_ds():
    from pylab import plot, legend, show, hold

    t=np.arange(-1,1,.01)[:,np.newaxis]
    x=np.sin(2*np.pi*t*.5+2)
    #xn=x + sin(2*pi*t*10)*.1
    xn=x+np.random.randn(*t.shape)*0.05
    out1 = np.zeros_like(x)[::10]
    out2 = np.zeros_like(x)[::10]

    
    downsample_1way(x,out1,ratio=10)
    downsample_2way(x,out2,ratio=10)

    plot(t,x,'c')
    hold(True)
    plot(t,xn,'k')
    plot(t[::10],out1,'r')
    plot(t[::10],out2,'g')

    legend(('original','noisy signal','1-way','2-way'))
    print "Notice phase lag from 1-way (causal) filter"
    hold(False)
    show()
    

if __name__ == '__main__':
    usage = """
    Downsample your data. Default ratio is 16, so N_out = floor(N_in/16).
    By default this script does a two-way filter (like filtfilt). To see why, try %prog -t
    If output file is not specified, it will have extension eeg (foo.dat -> foo.eeg)
    
    %prog your_dat_file.dat [output filename]
    %prog -h displays help"""
    parser = OptionParser(usage)
    parser.add_option("-1",help="One-way filter.",action="store_false",default=True,dest="two_way")
    parser.add_option("-r",help="Downsampling ratio",type=int,action="store",dest="ratio",default=RATIO)
    parser.add_option("-t",help="Test 1-way and two-way filtering.",action="store_true",dest="test",default=False)
    parser.add_option("--in-type",help="Datatype of input file. One of i2,i4,i8,f4,f8. Default i2 (16-bit integer).",action="store",dest="dt_in",default=DT_IN)
    parser.add_option("--out-type",help="Datatype of output file. One of i2,i4,i8,f4,f8. Default i2 (16-bit integer).",action="store",dest="dt_out",default=DT_OUT)
    (opts,args) = parser.parse_args()
    if opts.test:
        test_ds()
        exit()
    OUTNAME = args[1] if len(args) > 1 else switch_ext(args[0],"eeg")
    DT_IN = opts.dt_in
    DT_OUT = opts.dt_out
    downsample(args[0],OUTNAME,two_way=opts.two_way,ratio=opts.ratio)