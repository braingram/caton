import os,sys
import numpy as np

if not os.path.exists('/home/joschu/Data/fake32'): os.mkdir('/home/joschu/Data/fake32')
os.chdir('/home/joschu/Data/fake32')

DatPath = '/home/joschu/Data/d11221/d11221.002.dat'
n_ch_dat = 6
n_s = os.path.getsize(DatPath)//n_ch_dat//2

extra_chans = [0,3,1,2] # I ran into some bad luck: the intra unit happens to have spikes that
# are only on channels 0 and 3--the pair of channels that are not connected in the fake probe!
# So I reorder them.

X_sc = np.memmap(DatPath,dtype=np.int16,shape = (n_s,n_ch_dat),mode='r')
Y_sc = np.memmap('fake32.dat',dtype=np.int16,shape = (n_s-70000,40),mode='w+')

for ch_start,s_start,ch_intra in zip(xrange(0,32,4),xrange(0,80000,10000),xrange(32,40)):
    ch_stop = ch_start + 4
    s_stop = s_start-70000 or None
    Y_sc[:,ch_start:ch_stop] = X_sc[s_start:s_stop,extra_chans]
    Y_sc[:,ch_intra] = X_sc[s_start:s_stop,4]

del X_sc
del Y_sc
