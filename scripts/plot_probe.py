#!/usr/bin/env python

from caton.probe_stuff import plot_probe
from optparse import OptionParser
import caton.myopts

usage = """%prog your_probe_file.probe
%prog -h displays help"""

parser = OptionParser(usage)
parser.add_option(caton.myopts.output)
(opts,args) = parser.parse_args()
plot_probe(args[0],opts.output)

