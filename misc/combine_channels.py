#!/usr/bin/env python

import os,re,sys,numpy as np
from scanf import sscanf
from optparse import OptionParser
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from caton.utils_misc import indir,switch_ext
from caton.output import write_xml
join = os.path.join

usage = """%prog data_dir format_string output_file
%prog -h displays help

example:
%prog ~/Data/blah Jul8aCh{channel}_t{trial}.dat ~/blah/outfile.dat
"""

parser = OptionParser(usage)
parser.add_option("-d",action="store_true",dest="dry_run",default=False,help="dry run")
(opts,args) = parser.parse_args()
print args

data_dir = args[0]
pretty_fmt_str = args[1]
datfilename = os.path.abspath(args[2])
dry_run = opts.dry_run

format_str = re.sub("{channel}|{trial}","%d",pretty_fmt_str)
print format_str
in_dtype = np.dtype(">i2")
out_dtype = np.int16
in_data_dir = lambda filename: join(data_dir,filename)
n_bytes = in_dtype.itemsize
get_first = lambda filename: sscanf(filename,format_str)[0]
get_second = lambda filename: sscanf(filename,format_str)[1]
get_channel0,get_trial0 = (get_first,get_second) if pretty_fmt_str.index("channel") < pretty_fmt_str.index("trial") else (get_second,get_first)

def is_right_type(filename):
    try: 
        sscanf(filename,format_str)
        return True
    except Exception: 
        return False    

orig_dir = os.path.abspath(os.curdir)
os.chdir(data_dir)

all_files = os.listdir(data_dir)

data_files = filter(is_right_type,all_files)
data_files = sorted(data_files,key = lambda filename: (get_trial0(filename),get_channel0(filename)))    

channels = set(map(get_channel0,data_files))
trials = set(map(get_trial0,data_files))
n_channels = len(channels)
n_trials = len(trials)

channels_start_at = min(channels)
trials_start_at = min(trials)

if channels_start_at == 1:
    print "Channels start at 1. Subtracting 1 from every channel index."
    get_channel = lambda s: get_channel0(s)-1
else: get_channel = get_channel0

if trials_start_at == 1:
    print "Trials start at 1. Subtracting 1 from every trial index."
    get_trial = lambda s: get_trial0(s)-1
else: get_trial = get_trial0
    
n_channels = max(map(get_channel,data_files))+1 
n_trials = max(map(get_trial,data_files))+1 




print "%i channels, %i trials"%(n_channels, n_trials)

trial_lengths = np.zeros(n_trials,dtype=np.int32)
print "computing trial lengths..."
for filename in data_files:
    channel,trial = get_channel(filename),get_trial(filename)
    if channel == 1:
        trial_lengths[trial] = os.path.getsize(filename)//2
        print "trial length %i: %d"%(trial,trial_lengths[trial])

trial_starts = np.append(0,np.cumsum(trial_lengths))
n_samples_total = trial_starts[-1]

if not dry_run: Data_sc = np.memmap(datfilename,shape=(n_samples_total,n_channels),dtype=out_dtype,mode="w+")

for filename in data_files:
    if not dry_run: X_c = np.fromfile(filename,dtype=in_dtype)
    channel,trial = get_channel(filename),get_trial(filename)
    print "inserting data from %s"%filename
    if not dry_run: Data_sc[trial_starts[trial]:trial_starts[trial+1],channel] = X_c
    
if not dry_run: 
    del Data_sc
xmlpath = switch_ext(datfilename,"xml")
print "writing %s"%xmlpath
write_xml(n_channels,0,0,24000.,xmlpath)

os.chdir(orig_dir)