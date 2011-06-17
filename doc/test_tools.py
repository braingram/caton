import os,sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from caton.validate import print_isect_table,print_report_from_dirs
from caton.core import classify_from_raw_data,extract_intra_spikes,processed_basename,intra_basename
from caton.utils_misc import *
import caton.utils_misc

DataDir = "/home/joschu/Data/d11221"

###fast mode:
#import caton.core
#caton.core.extract_wave = caton.core.extract_wave_simple
#import caton.subset_sorting
#caton.subset_sorting.MINCLUSTERS = 5
#caton.subset_sorting.MAXCLUSTERS = 5

def test_and_write_report(DatFileName,ProbeFileName,IntraChannel,ResultDir,ReportFile,DontSort=False):
        print("Processing dat file %s with probe file %s"%(DatFileName,ProbeFileName))
        if not DontSort: classify_from_raw_data("batch",join(DataDir,DatFileName),ProbeFileName = ProbeFileName,output_dir=ResultDir)
        IntraDir = join(ResultDir,intra_basename(DatFileName))
        if os.path.exists(IntraDir):
                os.system("rm -rf %s"%IntraDir)
        extract_intra_spikes(join(DataDir,DatFileName),IntraChannel = IntraChannel,output_dir=ResultDir)
    
        caton.utils_misc.OUTPUT_FILE = ReportFile
        print_report_from_dirs(join(ResultDir,processed_basename(DatFileName,ProbeFileName,"batch")),IntraDir)
        oprint("")
        caton.utils_misc.OUTPUT_FILE = None
        
