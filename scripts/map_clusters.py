#!/usr/bin/env python

from optparse import OptionParser
import tables,os,shutil
from caton.utils_misc import *
from caton.output import write_clu
usage = """%prog blah.h5 "1 2 3" "4 4 5" """

parser = OptionParser(usage)
parser.add_option("--newdir",action="store",dest="newdir",default=False,
                help="""Store in new directory? Otherwise, modifies in place.""")

(opts,args) = parser.parse_args()

h5name = args[0]
basename = basename_noext(h5name)
newdir = opts.newdir
olddir = dirname(h5name)

oldchans = np.array([int(i) for i in args[1].split()],dtype=np.int32)
newchans = np.array([int(i) for i in args[2].split()],dtype=np.int32)

if newdir:
    shutil.copytree(dirname(abspath(h5name)),newdir)
    newfile = tables.openFile(join(newdir,h5name),"r+")
    newtable = newfile.root.SpikeTable
    oldtable = tables.openFile(h5name,"r").root.SpikeTable
else:
    newfile = tables.openFile(h5name,"r+")
    newtable = oldtable =   newfile.root.SpikeTable


n_clu = oldtable.cols.clu[:].max()+1
old2new = np.arange(n_clu,dtype=np.int32)
old2new[oldchans] = newchans

newtable.cols.clu[:] = old2new[oldtable.cols.clu[:]]
write_clu(newtable.cols.clu[:],switch_ext(join(newdir or olddir,h5name),"clu.1"))
