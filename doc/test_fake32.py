import os,sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from caton.utils_misc import *
from test_tools import test_and_write_report
join = os.path.join

DatFileName = "/home/joschu/Data/fake32/fake32.dat"
ResultDir = join("generated_data","fake32")
AuxDir = "aux"
mkdir_maybe_rm(ResultDir)

with open(join(ResultDir,"report.txt"),"w") as report:
    for IntraChannel in xrange(32,40):
        test_and_write_report(DatFileName,ProbeFileName=join(AuxDir,"fake32.probe"),ResultDir=ResultDir,ReportFile=report,IntraChannel=IntraChannel,
                              DontSort = (IntraChannel!=32))