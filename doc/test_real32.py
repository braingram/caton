import os,sys,re
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from caton.utils_misc import *
join = os.path.join

from time import time

DatFileName = "/home/joschu/Data/marn2dec13/marn2dec13.002-5/marn2dec13.002-5.dat"
ResultDir = join("generated_data","real32")

mkdir_maybe_rm(ResultDir)

batch_start = time()
os.system("../cluster_from_raw_data.py %s -n 200000 -o %s"%(DatFileName,ResultDir))
batch_end = time()

generalize_start = time()
BatchDir = join(ResultDir,filter(lambda string: string.endswith("batch"),os.listdir(ResultDir))[0])
os.system("../generalize_from_raw_data.py %s %s -o %s"%(DatFileName,BatchDir,ResultDir))
generalize_end = time()

with open(join(ResultDir,"times.txt"),"w") as fd:
    fd.write("""
    batch time: %i
    generalize time: %i
    """%(batch_end-batch_start,generalize_end-generalize_start))