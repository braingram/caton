from __future__ import division
import numpy as np
from core import *
import scipy.stats as ss

n_channels = 4
sample_rate = 25000.
time_length = 300
n_units = 5
f_min = 1; f_max = 4
T_min = .0004; T_max = .0006
a_min = 30.; a_max = 300.
noise_amp = 10.

data = np.zeros((int(time_length*sample_rate),n_channels),dtype=np.int32)
tseries = np.arange(0,time_length,1/sample_rate)    


f_arr = ss.uniform(f_min,f_max).rvs(n_units)
T_arr = ss.uniform(T_min,T_max).rvs(n_units)
a_arr = ss.uniform(a_min,a_max).rvs(n_units)
ch_arr = ss.randint(0,n_channels-1).rvs(n_units)    
n_spikes_arr = [ss.poisson(time_length*f).rvs() for f in f_arr]
spike_times = [ np.sort(ss.uniform(0,time_length).rvs(n_spikes)) for n_spikes in n_spikes_arr]

def add_spike(t,t_width,amp,ch):
    s = np.searchsorted(tseries,t)
    s_width = int(t_width*sample_rate)
    data[(s-s_width):(s+s_width),:] += (.75*amp*np.sin(-np.pi*(t-tseries[(s-s_width):(s+s_width)])/t_width)).reshape(-1,1)
    data[(s-s_width):(s+s_width),ch] += .25*amp*np.sin(-np.pi*(t-tseries[(s-s_width):(s+s_width)])/t_width)
    

for i_unit in xrange(n_units):
    for time in spike_times[i_unit]:
        add_spike(time,T_arr[i_unit],a_arr[i_unit],ch_arr[i_unit])            

white_noise = ss.uniform(-noise_amp,noise_amp).rvs(data.shape)
data += white_noise
       
fakedata_dir = os.path.join(caton.OUT_DIR,'fakedata')
if os.path.exists(fakedata_dir):
    os.system('rm -r '+fakedata_dir)
os.mkdir(fakedata_dir)
os.chdir(fakedata_dir)

import output
np.int16(data).tofile('fakedata1.dat')
write_info_about_file('fakedata1.dat',n_channels,'int16',sample_rate=sample_rate)
all_times = np.concatenate([np.int32(time_arr*sample_rate) for time_arr in spike_times])
all_clus = np.concatenate([[i]*len(spike_times[i]) for i in xrange(len(spike_times))])
sort_inds = np.argsort(all_times)
output.write_res(all_times[sort_inds],'fakedata1','')
output.write_clu(all_clus[sort_inds],'fakedata1','')