#!/usr/bin/env python

import sys,os

from optparse import OptionParser
from caton.core import extract_intra_spikes
import caton.myopts

usage = """%prog your_dat_file.dat intra_channel
%prog -h displays help"""

parser = OptionParser(usage)
parser.add_option(caton.myopts.output)
parser.add_option("-e",action="store",dest="ExtraStr",
                  help="""Optionally specify a number of extracellularly recorded channels to be included in spk file in the form
                  -e "0 1 2 3".""")



(opts,args) = parser.parse_args()

if len(args) != 2:
    parser.error("Wrong number of arguments")

DatFileName,IntraChannel = args
if not os.path.exists(DatFileName):
    parser.error("File not found: %s"%DatFileName)


ExtraChannels = [int(char) for char in opts.ExtraStr.split()] if opts.ExtraStr is not None else None
    
if __name__ == '__main__':
    extract_intra_spikes(DatFileName,IntraChannel,output_dir=opts.output,ExtraChannels = ExtraChannels)
