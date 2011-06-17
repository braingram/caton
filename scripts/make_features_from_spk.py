#!/usr/bin/env python

from caton.features import make_features_from_spk
from optparse import OptionParser

usage = """%prog your_spk_file.spk
%prog -h displays help"""

parser = OptionParser(usage)
(opts,args) = parser.parse_args()
make_features_from_spk(args[0])

