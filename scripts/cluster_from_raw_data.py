#!/usr/bin/env python
import sys,os,shutil
from optparse import OptionParser
from caton.core import classify_from_raw_data
import caton.myopts

usage = """
This is the main script that you use to spike-sort your data.
Just make the probe file and you're good to go (see documentation).
If you don't specify a probe file or probe directory, I will use every probe file in directory of dat file.
I will prompt you for sample rate and number of channels if no xml file is found.

%prog your_dat_file.dat [options]
%prog -h displays help"""
parser = OptionParser(usage)
parser.add_options([caton.myopts.probe,caton.myopts.max_spikes,caton.myopts.output,caton.myopts.fast,caton.myopts.fast2,caton.myopts.params])
parser.add_option("--probe-dir",action="store",dest="probe_dir",help="Run the clustering script once for every probe file in the specified directory. Defaults to the directory of the data file. If you want to select a single probe file, use the -p (or --probe=) option.")

if __name__ == '__main__':
    (opts,args) = parser.parse_args()                
        
    if len(args) == 0:
        parser.error("Must specify a dat file")
                        
    DatFileName = args[0]
    if not os.path.exists(DatFileName):
        parser.error("raw data file not found: %s"%DatFileName)        

    if opts.probe is not None:
        probe_files = [opts.probe]
    else:
        probe_dir = opts.probe_dir or os.path.dirname(os.path.abspath(DatFileName))
        print "I'll use all probe files in %s"%probe_dir
        probe_files = [os.path.join(probe_dir,fname) for fname in os.listdir(probe_dir) if fname.endswith(".probe")]
    print "%s will be run for each probe file in %s"%(__file__,probe_files)
    
    if len(probe_files) == 0: parser.error("no probe files found!")
    for probe_file in probe_files:
        classify_from_raw_data("batch",DatFileName,probe_file,max_spikes=opts.max_spikes,output_dir=opts.output)
