#!/usr/bin/env python

from caton.features import plot_features
from optparse import OptionParser
import caton.myopts

usage = """%prog number_of_features"""

parser = OptionParser(usage)
parser.add_option(caton.myopts.output)
parser.add_option("-n",action="store",dest="F",default=3,
                help="""Number of features.""")
parser.add_option("-s","--savefig",action="store_true",dest="savefig",default=False,
                help="""Save figure as .png""")

(opts,args) = parser.parse_args()
plot_features(opts.F,output_dir=opts.output,savefig=opts.savefig)
