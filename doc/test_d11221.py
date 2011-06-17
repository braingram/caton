import os,sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from caton.utils_misc import *
from test_tools import test_and_write_report

if not os.path.exists("generated_data"): os.mkdir(generated_data)
DataDir = "/home/joschu/Data/d11221"
ResultDir = join("generated_data","d11221")
AuxDir = "aux"
mkdir_maybe_rm(ResultDir)

# filenames of probe files
just_tet_setup = {"ProbeFileName":join(AuxDir,"just_tet.probe"),"IntraChannel":4}
fil_tet_setup = {"ProbeFileName":join(AuxDir,"fil_tet.probe"),"IntraChannel":5}
tet_fil_fil_setup = {"ProbeFileName":join(AuxDir,"tet_fil_fil.probe"),"IntraChannel":6}

DatInfo = {
    "d11221.001.dat":just_tet_setup,
    "d11221.002.dat":just_tet_setup,
    "d1122101.dat":tet_fil_fil_setup,
    "d1122102.dat":tet_fil_fil_setup,
    "d1122103.dat":tet_fil_fil_setup,    
    "d1122104.dat":fil_tet_setup,
    "d1122105.dat":fil_tet_setup,
    "d1122106.dat":fil_tet_setup,
    "d1122107.dat":fil_tet_setup,
    "d1122108.dat":tet_fil_fil_setup,
    "d1122109.dat":tet_fil_fil_setup}

with open(join(ResultDir,"report.txt"),"w") as report:
    for DatFileName,Info in DatInfo.items():
        test_and_write_report(DatFileName,ResultDir=ResultDir,ReportFile=report,**Info)
    
# The data files in d11221 have a tetrode, plus zero, one, or two filtered channels
# There are three different setups, which correspond to the three probe files in this directory

# d11221.001 tet 0-3 intra 4
# d11221.002 tet 0-3 intra 4
# d1122101 tet 0-3 fil 4-5 intra 6
# d1122102 tet 0-3 fil 4-5 intra 6
# d1122103 tet 0-3 fil 4-5 intra 6
# d1122104 fil 0 tet 1-4 intra 5
# d1122105 fil 0 tet 1-4 intra 5
# d1122106 fil 0 tet 1-4 intra 5
# d1122107 fil 0 tet 1-4 intra 5
# d1122108 tet 0-3 fil 4-5 intra 6
# d1122109 tet 0-3 fil 4-5 intra 6