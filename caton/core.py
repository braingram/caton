from __future__ import with_statement, division
import itertools as it, numpy as np, scipy.signal as signal
from scipy.stats import rv_discrete
from scipy.stats.mstats import mquantiles
from xml.etree.ElementTree import ElementTree
import re, tables, json, os

from tables import IsDescription,Int32Col,Float32Col,Int8Col
import probe_stuff, output
from output import get_pars_from_xml, get_pars_from_prompt, get_dat_pars
from utils_graphs import contig_segs,complete_if_none
from utils_misc import is_numerical, is_bool, find_file_with_ext, consumerize, indir, dump, to2d, tofloat32, basename_noext, memoized, get_padded, switch_ext, naive_maximize, inrange, dict_append, flatdist
from CEM_extensions import bincount, class_means, class_covs, class_wts, sqnorm, subset_inds, connected_components
from features import get_features,compute_pcs
from subset_sorting import cluster_withsubsets, spike_subsets
from progressbar import ProgressBar,SimpleProgress,Bar,Percentage
from interp_stuff import interp_around_peak
from os.path import join,abspath,dirname
from time import sleep

### Placeholders. These parameters are set in PARAMS.py
T_BEFORE = T_AFTER = T_JOIN_CC = F_LOW = THRESH_SD = DETECT_POSITIVE = DTYPE = SEPARATE_CHANNELS_PCA = REGET_FEATURES = SORT_CLUS_BY_CHANNEL = FPC = BUTTER_ORDER = CHUNK_SIZE = CHUNK_OVERLAP = CHUNKS_FOR_THRESH = INTERP_METHOD = None

execfile(join(dirname(abspath(__file__)),"PARAMS.py"))

### These parameters will be set by set_globals_samples
SAMPLE_RATE = F_HIGH = None # sample rate (Hz), top of passband (Hz)
S_BEFORE = S_AFTER = S_TOTAL = S_JOIN_CC = None

PC_3s = None # feature vectors. set by load_features(). not necessarily 3 despite name.
PC_df = None # set by get_features_allch()
N_CH = None

PARAMETERS = dict([(k,v) for k,v in globals().items() if (is_numerical(v) or is_bool(v))])

def set_globals_samples(sample_rate):
    """parameters are set in terms of time (seconds). this sets corresponding parameters in terms of sample rate. should be run before any processing"""
    
    global SAMPLE_RATE, F_HIGH
    global S_BEFORE, S_AFTER, S_TOTAL, S_JOIN_CC

    SAMPLE_RATE = sample_rate
    F_HIGH = .95*SAMPLE_RATE/2
    S_BEFORE = int(T_BEFORE*SAMPLE_RATE)
    S_AFTER = int(T_AFTER*SAMPLE_RATE)
    S_TOTAL = S_BEFORE + S_AFTER
    S_JOIN_CC = T_JOIN_CC*SAMPLE_RATE
    
def print_params():
    print "PARAMETERS"
    print "----------"
    print "T_BEFORE: %s"%T_BEFORE
    print "T_AFTER: %s"%T_AFTER
    print "F_LOW (Hz): %s"%F_LOW
    print "THRESH_SD: %s"%THRESH_SD    
    
####################################
######## High-level scripts ########
####################################

def classify_from_raw_data(JobType,DatFileName,ProbeFileName,max_spikes=None,output_dir=None,clu_dir=None):
    """Top level function that starts a data processing job. JobType is "batch" or "generalize". """
    print_params()

    if not os.path.exists(DatFileName):
        raise Exception("Dat file %s does not exist"%DatFileName)
    DatFileName = os.path.abspath(DatFileName)
    
    n_ch_dat,sample_rate = get_dat_pars(DatFileName)
    set_globals_samples(sample_rate)
    
    DatDir = os.path.abspath(os.path.dirname(DatFileName))
    probe_stuff.load_probe(ProbeFileName)
    global N_CH; N_CH = probe_stuff.N_SITES
    
    basename = processed_basename(DatFileName,ProbeFileName,JobType)
    OutDir = join(output_dir,basename) if output_dir else join(DatDir,basename)
    if clu_dir is not None: clu_dir = os.path.abspath(clu_dir)
    with indir(OutDir):    
        Channels_dat = [site.dat for site in probe_stuff.PROBE_SITES]
        if JobType == "batch":
            cluster_from_raw_data(basename,DatFileName,n_ch_dat,Channels_dat,probe_stuff.PROBE_GRAPH,max_spikes)
        realOutDir = os.getcwd()
        #elif JobType == "generalize":
            #generalize_group_from_raw_data_splitprobe(basename,DatFileName,n_ch_dat,Channels_dat,probe_stuff.PROBE_GRAPH,max_spikes,clu_dir)                        
    return realOutDir # return the directory where the clustering data was just written

#def spike_dtype():
  #  return np.dtype([("time",np.int32),("st",np.int8,N_CH),("wave",np.float32,(S_TOTAL,N_CH)),("fet",np.float32,(N_CH,FPC)),("clu",np.int32)])            

def spike_dtype():
    class description(IsDescription):    
        time = Int32Col()
        st = Int8Col(shape=(N_CH,))
        wave = Float32Col(shape=(S_TOTAL,N_CH))
        fet = Float32Col(shape=(N_CH,FPC))
        clu = Int32Col()
    return description

def make_table(filename,tablename):
    h5file = tables.openFile(filename,mode="w",title="Spike info file")
    new_table = h5file.createTable(h5file.root,name = tablename,description=spike_dtype())
    return new_table,h5file
    
    
            
def cluster_from_raw_data(basename,DatFileName,n_ch_dat,Channels_dat,ChannelGraph,max_spikes):
    """Filter, detect, extract, and cluster on raw data."""
    ### Detect spikes. For each detected spike, send it to spike writer, which writes it to a spk file.
    ### List of times is small (memorywise) so we just store the list and write it later.

    hdf5file = tables.openFile(basename+".h5","w")
    spike_table = hdf5file.createTable("/","SpikeTable_temp",spike_dtype())
    np.savetxt("dat_channels.txt",Channels_dat,fmt="%i")
    hdf5file.createArray("/","DatChannels",Channels_dat)
    
    for Spk,PeakSample,ST in extract_spikes(DatFileName,n_ch_dat,Channels_dat,ChannelGraph,max_spikes):
        spike_table.row["wave"] = Spk
        spike_table.row["time"] = PeakSample
        spike_table.row["st"] = ST
        spike_table.row.append()
        
    spike_table.flush()    

    
    ### Feature extraction on spikes    
    if SEPARATE_CHANNELS_PCA:
        if REGET_FEATURES: reget_features(spike_table.cols.wave[:10000])
        else: load_features()
        # work around pytables bug
        #spike_table.cols.fet[:] = [project_features(wave) for wave in spike_table.cols.wave]
        for row in spike_table: 
            row["fet"] = project_features(row["wave"])
            row.update()
    else:
        get_features_allch(spike_table.cols.wave[:10000])
        spike_table.cols.fet[:] = [project_features_allch(wave) for wave in spike_table.cols.wave]        
                    
    ### And use batch clustering procedure (EM) to get clusters.
    CluArr = cluster_withsubsets(spike_table,SORT_CLUS_BY_CHANNEL)
    spike_table.cols.clu[:] = CluArr
    TmArr = spike_table.cols.time[:]

    final_table = hdf5file.createTable("/","SpikeTable",spike_dtype())
    for ind in good_inds(CluArr,TmArr):
        final_table.append(spike_table[ind:(ind+1)])
    #final_table.append(spike_table[good_inds(CluArr,TmArr)])
    final_table.flush()

    hdf5file.removeNode("/SpikeTable_temp")
    
    klusters_files(final_table,basename)
    hdf5file.close()
    

def good_inds(CluArr,TmArr):
    """Returns indices of all spikes that (1) do not appear to be a copy of another spike, and (2) are not noise (i.e., cluster 0)"""
    return [i_row for i_row in TmArr[:-1].argsort() if 
            CluArr[i_row] != 0 # not in cluster with discarded spikes
            and not (
                (TmArr[i_row+1]-TmArr[i_row] < S_JOIN_CC) and # same time
                (CluArr[i_row+1] == CluArr[i_row]))] # same cluster
    
#def generalize_group_from_raw_data_splitprobe(basename,DatFileName,n_ch_dat,Channels_dat,ChannelGraph,max_spikes,clu_dir):
    #"""Filter, detect, extract, using fet and clu files in clu_dir to generate model."""    
    
    #SpkFileName = basename+'.spk.1'
    #OldDir = clu_dir
    #OldCluList = output.read_clu(find_file_with_ext(OldDir,'.clu.1'))
    #OldFetList = output.read_fet(find_file_with_ext(OldDir,'.fet.1'))[:,:(N_CH*FPC)]
    #STList = np.load(join(OldDir,"ST.npy"))
    
    #### Detect spikes. For each detected spike, send it to spike writer, which writes it to a spk file.
    #### List of times is small (memorywise) so we just store the list and write it later.
    #load_features()
    #SpkWriter = consumerize(spk_writer,"SpkStream",SpkFileName = SpkFileName)
    #SpkSorter = consumerize(generalize_spk_cuts,"NewSpkStream","NewSTStream",
                            #OldFetList = OldFetList,OldCluList = OldCluList,OldSTList = STList,ChannelGraph=ChannelGraph)
    #TmList,CluList,FetList,STList=[],[],[],[]
    #for Spk,PeakSample,ST in extract_spikes(DatFileName,n_ch_dat,Channels_dat,ChannelGraph,max_spikes):
        #Clu,Spk,Fet = SpkSorter.send(Spk,ST)
        #CluList.append(Clu)
        #SpkWriter.send(Spk)
        #FetList.append(Fet)
        #TmList.append(PeakSample)        
        #STList.append(ST)
    #SpkWriter.close()
    #SpkSorter.close()

    ### Write all the files we need for klusters.
    #write_files(basename,CluList,TmList,to2d(np.array(FetList)),ChannelGraph)
    
def combine_h5s(*dirs):
    "combine the data from a bunch of h5 files. Also make klusters files"    
    outdirname = common_start(dirs) + "COMBINED"
    print dirs
    print os.path.abspath(os.curdir)
    h5files = [tables.openFile(find_file_with_ext(dir,".h5"),mode="r") for dir in dirs]
    spike_tables = [h5file.root.SpikeTable for h5file in h5files]
    with indir(outdirname):
                
        global N_CH,S_TOTAL,FPC,SAMPLE_RATE
        num_ch =[table.cols.st.shape[1] for table in spike_tables]
        N_CH = max(num_ch)
        S_TOTAL = table.cols.wave.shape[1]
        FPC = table.cols.fet.shape[2]
        new_file = tables.openFile(outdirname+".h5",mode="w")
        new_table = new_file.createTable("/","SpikeTable",spike_dtype())
        clu_start = np.arange(0,5000,100)
        SAMPLE_RATE = 25000. # Doesn't actually matter for this script, but write_xml wants it
       
        # files in same order as tables
        clu2DatChannels = {}
        for (i_file,h5file) in enumerate(h5files):
            for clu in xrange(clu_start[i_file],clu_start[i_file+1]):
                clu2DatChannels[clu] = list(h5file.root.DatChannels) 
        dump("clu2DatChannels.pickle",clu2DatChannels)                        
        
        triples = [(row["time"],i_spike,i_table) for (i_table,table) in enumerate(spike_tables) for (i_spike,row) in enumerate(table)]
        for time,i_spike,i_table in sorted(triples, key = lambda tup: tup[0]):
            oldrow = spike_tables[i_table][i_spike]
            new_table.row["time"] = time
            new_table.row["fet"]= zero_pad(oldrow["fet"],N_CH)
            new_table.row["st"] = zero_pad(oldrow["st"],N_CH)
            new_table.row["wave"] = zero_pad(oldrow["wave"],N_CH)
            new_table.row["clu"] = oldrow["clu"]+clu_start[i_table]
            new_table.row.append()
        new_table.flush()
        klusters_files(new_table,outdirname)
                
        new_file.close()

def zero_pad(arr,newlen):
    n = arr.shape[0]
    if n == newlen: return arr
    elif n < newlen: return np.concatenate(  (arr,np.zeros((newlen-n,)+arr.shape[1:],dtype=arr.dtype))  ,0)
        
        
        
def common_start(string_list):
    for n in range(len(string_list[0])):
        substring = string_list[0][:n]
        if not all((substring == string[:n] for string in string_list)): break
    return substring[:-1]

        
def nonduplicates(TmArr, CluArr):
    return np.flatnonzero((TmArr[1:]-TmArr[:-1] > S_JOIN_CC) | (CluArr[1:]!=CluArr[:-1]))

def duplicates(TmArr, CluArr):
    return np.flatnonzero((TmArr[1:]-TmArr[:-1] < S_JOIN_CC) & (CluArr[1:]==CluArr[:-1]))


def klusters_files(table,basename):
    CluFileName,FetFileName,ResFileName,SpkFileName,XMLFileName = (basename + ext for ext in [".clu.1",".fet.1",".res.1",".spk.1",".xml"])
    output.write_clu(table.cols.clu[:],CluFileName)
    output.write_fet(to2d(table.cols.fet[:]),FetFileName,samples=table.cols.time[:])
    output.write_res(table.cols.time[:],ResFileName)
    output.write_spk(table.cols.wave[:],SpkFileName)
    output.write_xml(n_ch=N_CH,n_samp=S_TOTAL,n_feat=N_CH*FPC,sample_rate=SAMPLE_RATE,filepath=XMLFileName)    


def write_files(basename,CluList,TmList,FetList,ChannelGraph=None,STArr=None):
    """Writes files that result from a clustering job."""
    CluFileName = basename+'.clu.1'
    FetFileName = basename+'.fet.1'
    ResFileName = basename+'.res.1'
    XMLFileName = basename+'.xml'    
    if CluList is not None: output.write_clu(np.array(CluList),CluFileName)
    output.write_res(np.array(TmList),ResFileName)
    output.write_fet(np.array(FetList),FetFileName,samples=np.array(TmList))
    output.write_xml(n_ch=N_CH,n_samp=S_TOTAL,n_feat=N_CH*FPC,sample_rate=SAMPLE_RATE,filepath=XMLFileName)
                
    if STArr is not None:
        np.save("ST.npy",STArr)
        
    #with open("params.json","w") as fd: json.dump(dict_append(PARAMETERS,dict(SAMPLE_RATE=str(SAMPLE_RATE),N_CH=N_CH)),fd)



def extract_intra_spikes(DatFileName,IntraChannel,output_dir=None,ExtraChannels=None):
    """extracts spikes times from intracellular data"""
    
    THRESH_FRAC = .5
    
    DatFileName = os.path.abspath(DatFileName)
    DatDir = os.path.dirname(DatFileName)
    basename = intra_basename(DatFileName)
    
    OutDir = join(output_dir,basename) if output_dir else join(DatDir,basename)
    
    with indir(OutDir):
        SpkFileName = basename+'.spk.1'
        
        n_ch_dat,sample_rate = get_dat_pars(DatFileName)
        global N_CH
        N_CH = 1 if ExtraChannels is None else 1 + len(ExtraChannels)
        set_globals_samples(sample_rate)
    
        print("extracting intracellular spikes from %s"%DatFileName)
    
        n_samples = num_samples(DatFileName,n_ch_dat,n_bytes=np.nbytes[DTYPE])
        AllDataArr = np.memmap(DatFileName,dtype=np.int16,shape=(n_samples,n_ch_dat),mode='r')
    
        b,a = signal.butter(3,100./(SAMPLE_RATE/2),'high') #filter at 100 Hz
        IntraArr = AllDataArr[:,IntraChannel].copy()
        IntraArr = signal.filtfilt(b,a,IntraArr)
        Thresh = IntraArr.max()*THRESH_FRAC    
    
        Segs = contig_segs(np.flatnonzero(IntraArr > Thresh),padding=2)
        TmList = map(lambda Seg: Seg[IntraArr[Seg].argmax()],Segs)
        CluList = np.ones(len(TmList),dtype=np.int)
        FetList = np.zeros((len(TmList),1),dtype=np.int)
        SpkList = [get_padded(IntraArr,PeakSample-S_BEFORE,PeakSample+S_AFTER) for PeakSample in TmList]
        SpkArr = np.array(SpkList)[:,:,np.newaxis]
        
        
        if ExtraChannels is not None:
            ExtraArr = AllDataArr[:,ExtraChannels].copy()
            #b,a = signal.butter(BUTTER_ORDER,(F_LOW/(SAMPLE_RATE/2),.95),'pass')                
            ExtraArr = filtfilt2d(b,a,ExtraArr)
            ExtraSpkList = [get_padded(ExtraArr,PeakSample-S_BEFORE,PeakSample+S_AFTER) for PeakSample in TmList]        
            ExtraSpkArr = np.array(ExtraSpkList)
            SpkArr *= ExtraSpkArr[0].max()/SpkArr[0].max()                       
            SpkArr = np.concatenate((np.array(ExtraSpkList),SpkArr),axis=2)
    
        output.write_spk(np.array(SpkArr),SpkFileName)
        write_files(basename,CluList,TmList,FetList,[])
    

    
###########################################################
############# Spike extraction helper functions ###########    
###########################################################


def extract_spikes(DatFileName,n_ch_dat,ChannelsToUse,ChannelGraph,max_spikes=None):
    
    n_samples = num_samples(DatFileName,n_ch_dat)
    b,a = signal.butter(BUTTER_ORDER,(F_LOW/(SAMPLE_RATE/2),F_HIGH/(SAMPLE_RATE/2)),'pass')    
    ProgBar = spikes_and_samples_bar(max_spikes,n_samples)
    
    with open(DatFileName,'r') as fd:
        # Use 5 chunks to figure out threshold
        n_samps_thresh = min(CHUNK_SIZE*CHUNKS_FOR_THRESH,n_samples)
        DatChunk = np.fromfile(fd,dtype=DTYPE,count=n_samps_thresh*n_ch_dat).reshape(n_samps_thresh,n_ch_dat)[:,ChannelsToUse]
        fd.seek(0)
        FilteredChunk = filtfilt2d(b,a,DatChunk.astype(np.int32))    
        Threshold = THRESH_SD*np.median(np.abs(FilteredChunk),axis=0)/.6745    
        
        spike_count = 0
        for s_start,s_end,keep_start,keep_end in chunk_bounds(n_samples,CHUNK_SIZE,CHUNK_OVERLAP):
            DatChunk = np.fromfile(fd,dtype=DTYPE,count=(s_end-s_start)*n_ch_dat).reshape(s_end-s_start,n_ch_dat)[:,ChannelsToUse]
            fd.seek(fd.tell()-CHUNK_OVERLAP*n_ch_dat*np.nbytes[DTYPE])
            FilteredChunk = filtfilt2d(b,a,DatChunk.astype(np.float32))    
            BinaryChunk = (FilteredChunk < -Threshold).astype(np.int8) if DETECT_POSITIVE else np.abs(FilteredChunk < -Threshold).astype(np.int8)
            IndListsChunk = connected_components(BinaryChunk,complete_if_none(ChannelGraph,N_CH),S_JOIN_CC)              
            for IndList in IndListsChunk:
                wave,s_peak,st = extract_wave(IndList,FilteredChunk)
                if inrange(s_start+s_peak,keep_start,keep_end):
                    spike_count += 1
                    yield wave,s_start + s_peak,st
            ProgBar.update(spike_count,s_end)
            if max_spikes is not None and spike_count >= max_spikes: break
        ProgBar.finish()
         
    
def chunk_bounds(n_samples,chunk_size,overlap):
    s_start = 0
    s_end = chunk_size
    keep_start = s_start
    keep_end = s_end-overlap//2
    yield s_start,s_end,keep_start,keep_end
    
    while s_end-overlap+chunk_size < n_samples:
        s_start = s_end-overlap
        s_end = s_start+chunk_size
        keep_start = keep_end
        keep_end = s_end-overlap//2
        yield s_start,s_end,keep_start,keep_end
        
    s_start = s_end-overlap
    s_end = n_samples
    keep_start = keep_end
    keep_end = s_end
    yield s_start,s_end,keep_start,keep_end   


#def generalize_spk_cuts(OldFetList,OldCluList,OldSTList,NewSpkStream,NewSTStream,ChannelGraph):

    #Mean_mf,Weight_m,Cov_mff = get_cluster_params(np.array(OldFetList)[:,:(N_CH*FPC)],OldCluList,ChannelGraph)
    #Cutoff_m = find_cutoffs(reshape_sep_chans(np.array(OldFetList)),np.array(OldCluList),np.array(OldSTList))
    #print("Cutoffs: %s"%Cutoff_m)
    
    #Fet_nc3 = it.imap(project_features,NewSpkStream)
    #global ChSubsets # global because subset_candidate_inds is memoized and just uses ST_c
    #ChSubsets = probe_stuff.ch_subsets()    
    #GoodCluClassifier = consumerize(best_good_clus,"Fet_nc3","ST_nc",Mean_mf=Mean_mf,Weight_m=Weight_m,Cov_mff = Cov_mff,ChSubsets=ChSubsets)
    #for Spk,ST in it.izip(NewSpkStream,NewSTStream):
        #fet_c3 = project_features(Spk)
        #Clu,LogL = GoodCluClassifier.send(fet_c3,ST)
        #if LogL < Cutoff_m[Clu]: Clu = 0            
        #yield Clu,Spk,fet_c3

def extract_wave_simple(IndList,FilteredArr):    
    IndArr = np.array(IndList,dtype=np.int32)
    SampArr = IndArr[:,0]
    ChArr = IndArr[:,1]
    
    PeakSample = SampArr[FilteredArr[SampArr,ChArr].argmin()]
    Wave = get_padded(FilteredArr,PeakSample-S_BEFORE,PeakSample+S_AFTER)
    return Wave,PeakSample,bincount(ChArr,N_CH).astype(np.bool8)

def extract_wave_interp(IndList,FilteredArr):
    IndArr = np.array(IndList,dtype=np.int32)
    SampArr = IndArr[:,0]
    ChArr = IndArr[:,1]
    
    PeakInd = FilteredArr[SampArr,ChArr].argmin()
    PeakSample,PeakChannel = SampArr[PeakInd],ChArr[PeakInd]
    WavePlus = get_padded(FilteredArr,PeakSample-S_BEFORE-1,PeakSample+S_AFTER+1)
    Wave = interp_around_peak(WavePlus,S_BEFORE+1,PeakChannel,S_BEFORE,S_AFTER,kind=INTERP_METHOD)                
    
    return Wave,PeakSample,bincount(ChArr,N_CH).astype(np.bool8)
    
extract_wave = extract_wave_interp

def filtfilt2d(b,a,x):
    out_arr = np.zeros_like(x)
    for i_ch in range(x.shape[1]):
        out_arr[:,i_ch] = signal.filtfilt(b,a,x[:,i_ch]) 
    return out_arr    
    
##################################################
####### Percentile stuff for generalize ##########
##################################################

#def find_cutoffs(Fet_nc3,Clu_n,ST_nc):
    #"""Return the likelihood cutoff for each cluster. Each incoming point is matched to
    #the most likely cluster. If the likelihood is lower than this cutoff, it is classified
    #as a noise point instead."""
    #global ChSubsets
    #ChSubsets = probe_stuff.ch_subsets()
    
    ## get cluster parameters needed for likelihood
    #GoodInds = np.flatnonzero(Clu_n != 0)
    #Mean_mf,Weight_m,Cov_mff = get_cluster_params(
        #to2d(Fet_nc3[GoodInds]),
        #Clu_n[GoodInds],
        #probe_stuff.PROBE_GRAPH)
    
    ## compute likelihood of every point in every cluster    
    #GoodClu_n,LogL_n = best_good_clus_batch(Fet_nc3,ST_nc,Mean_mf, Weight_m,Cov_mff,ChSubsets)    
    #M = Mean_mf.shape[0]
    #IndList_m = subset_inds(GoodClu_n,M=M)
    #Cutoff_m = np.zeros(M,dtype=np.float32)    
    #for m in xrange(1,M):
        #GoodInds_in_m = np.array(IndList_m[m],dtype=np.int32)[Clu_n[IndList_m[m]] != 0]
        #BadInds_in_m = np.array(IndList_m[m],dtype=np.int32)[Clu_n[IndList_m[m]] == 0]
        
        #if len(BadInds_in_m) > 0:
            #GoodDist = bootstrap_rv(LogL_n[GoodInds_in_m])
            #BadDist = bootstrap_rv(LogL_n[BadInds_in_m])
            #Cutoff_m[m] = naive_maximize(LossFunction,mquantiles(LogL_n[GoodInds_in_m],prob=np.arange(0,.25,.005)),func_args=dict(GoodDist=GoodDist,BadDist=BadDist),MinimizeInstead=True)                      
        #else:
            #Cutoff_m[m] = -np.inf
            
    #del ChSubsets
    #return Cutoff_m
    
#def LossFunction(cutoff,GoodDist,BadDist):
    #false_pos = 1-BadDist.cdf(cutoff)
    #false_neg = GoodDist.cdf(cutoff)
    #return false_pos+false_neg

#def bootstrap_rv(values):
    #return rv_discrete(name="bootstrap_dist",values=(values,flatdist(values.size)))
    
#def select_highamp_group(Fet_c3,ST_c):
    #Amp_c = Fet_c3[:,0]**2
    #SubsetCandInds = subset_cand_inds(ST_c)    
    #Subset_amp = [Amp_c[ChSubsets[i_subset]].sum() for i_subset in SubsetCandInds]
    #return SubsetCandInds[np.argmax(Subset_amp)]

#@memoized
#def subset_cand_inds(ST_c):
    #NumST = [ST_c[subset].sum() for subset in ChSubsets]
    #max_num = max(NumST)
    #return np.flatnonzero(NumST == max_num)

#def asarray(mat):
    #return mat.toarray() if mat.__module__ == 'scipy.sparse.csr' else mat.A

#def best_good_clus(Fet_nc3,ST_nc, Mean_mf, Weight_m, Cov_mff,ChSubsets):
    #C = N_CH; M = Mean_mf.shape[0]
    
    #Mean_nc3 = Mean_mf.reshape(M,C,FPC)
                
    #Cov_mc3c3 = [Cov_mff[m].reshape((C,FPC,C,FPC)) for m in xrange(M)]
    #Mean_slices = [to2d(Mean_nc3[:,channels,:]) for channels in ChSubsets]
    #InvSqrtCov_slices = [[np.asmatrix(np.linalg.inv(Cov_mc3c3[m][np.ix_(channels,range(FPC),channels,range(FPC))].reshape(FPC*len(channels),FPC*len(channels)))) for m in xrange(M)] for channels in ChSubsets]

    #for fet_c3,st_c in it.izip(Fet_nc3,ST_nc):
        #i_group = select_highamp_group(fet_c3,st_c)
        #channels = ChSubsets[i_group]
        #LogP_nm = compute_LogP(fet_c3[channels,:].reshape(1,-1),Mean_slices[i_group],Weight_m,InvSqrtCov_slices[i_group])
        #yield LogP_nm.argmax(),LogP_nm.max()
        
        
#def best_good_clus_batch(Fet_nc3,ST_nc,Mean_mf, Weight_m,InvSqrtCov_mff,ChSubsets):
    #"""
    #Returns
    #------
    #BestGoodClu_n: m that maximizes LogP_nm on selected channel subset
    #LogL_n: value of LogP at maximum"""
    #BestGoodClu_n,LogL_n = zip(*list(best_good_clus(Fet_nc3,ST_nc,Mean_mf, Weight_m,InvSqrtCov_mff,ChSubsets)))
    #return np.array(BestGoodClu_n),np.array(LogL_n)


###########################
############# I/O #########
###########################

#def infer_dtype(DatFileName):
    #dtypes = [np.dtype("<i2"),np.dtype(">i2")]
    #avg_mags = map(lambda dtype:
                  #np.log(1+np.abs(np.fromfile(DatFileName,dtype=dtype,count=100)).sum()),
                  #dtypes)
    #return dtypes[np.argmin(avg_mags)]


#def spk_writer(SpkStream,SpkFileName):
    #with open(SpkFileName,'w') as fd:        
        #while True:
            #spk = SpkStream.next()
            #spk.astype(np.int16).tofile(fd)
            #yield

        

#########################################
#### Feature Extraction #################
#########################################

def shift_waves(X_sc,IntShifts):
    return np.array([np.roll(X_sc,shift,axis=0) for shift in IntShifts])

def features_maxminsum(X_sc):
    return np.hstack((X_sc.max(axis=0),X_sc.min(axis=0),X_sc.sum(axis=0)))

def load_features():
    global PC_3s
    PC_3s = get_features(S_BEFORE,S_AFTER,SAMPLE_RATE,FPC)
    #import matplotlib.pyplot as plt
    #for i in xrange(3):
    #    plt.plot(PC_3s[i])
    #plt.show()

def get_features_allch(X_nsc):
    global PC_df
    PC_df = compute_pcs(to2d(X_nsc))[:(FPC*N_CH)]    
def reget_features(X_nsc):
    global PC_3s
    PC_3s = compute_pcs(X_nsc[:,:,0])[:FPC]
    import matplotlib.pyplot as plt
    for i in xrange(3):
        plt.plot(PC_3s[i])
    plt.show()


def project_features_allch(X_sc):
    return 100*np.dot(PC_df,X_sc.flatten())
def nproject_features_allch(X_nsc):
    return 100*np.tensordot(to2d(X_nsc),PC_df,axes = [1,1])

def project_features(X_sc):
    return 100.*np.dot(PC_3s,X_sc).T
def project_features_flat(X_sc):
    return project_features(X_sc).T.flatten()
def nproject_features(X_nsc):
    return 100.*np.tensordot(X_nsc,PC_3s,axes=[1,1])
def nproject_features_flat(X_nsc):
    Fet_nc3 = project_features(X_nsc)
    return to2d(Fet_nc3)


def pca_components(X_ns,F):
    print("solving for pca components")
    X_ns = X_ns.astype(np.float64)
    CentX_ns = X_ns - X_ns.mean(axis=0)
    CompVecs_ss = np.linalg.svd(CentX_ns)[2]
    return tofloat32(CompVecs_ss[:F])




#########################################
######## Misc helper functions ##########
#########################################
    
def processed_basename(DatFileName,ProbeFileName,JobType):
    return "%s_%s_%s"%(basename_noext(DatFileName),basename_noext(ProbeFileName),JobType)

def intra_basename(DatFileName):
    return "%s_%s"%(basename_noext(DatFileName),"intracellular")                
    
#def get_cluster_params(FetList,CluList,ChannelGraph):
    #FetArr = np.array(FetList,dtype=np.float32)
    #CluArr = np.array(CluList,dtype=np.int32)
    #N,F = FetArr.shape
    #M = CluArr.max()+1
    
    #Mean_mf, Cov_mff, Weight_m = m_step(CluArr,FetArr,F,M,N)
    #return Mean_mf,Weight_m,Cov_mff    

#def m_step(Class_n,X_nf,F,M,N):
    #Mean_mf = class_means(X_nf,Class_n,M)
    #Cov_mff = class_covs(X_nf,Mean_mf,Class_n,M)
    #Weight_m = class_wts(Class_n,M)
    #return Mean_mf, Cov_mff, Weight_m

#def compute_LogP(X_nf, Mean_mf, Weight_m, InvSqrtCov_mff):
    #N = X_nf.shape[0]
    #M,F = Mean_mf.shape
    
    #LogP = np.zeros((N,M),dtype=np.float32)
    #for m in xrange(M):
        #Vec2Mean_nf = Mean_mf[m] - X_nf
        #LogInvSqrtDet = np.log(InvSqrtCov_mff[m].diagonal()).sum()        
        #Mahal_n = sqnorm(InvSqrtCov_mff[m]*Vec2Mean_nf.T)
        #LogP[:,m] = - Mahal_n/2 + LogInvSqrtDet + np.log(Weight_m[m]) - np.log(2*np.pi)*F/2 ###+ logrootdet because we're using the inverse
    #return LogP     



#def reshape_sep_chans(Arr_ijkf):
    #return Arr_ijkf.reshape(Arr_ijkf.shape[:-1] + (Arr_ijkf.shape[-1]/FPC,FPC))
    
def num_samples(FileName,n_ch_dat,n_bytes=2):
    total_bytes = os.path.getsize(FileName)
    if total_bytes % (n_ch_dat*n_bytes) != 0:
        raise Exception("Size of file %s is not consistent with %i channels and %i bytes"%(FileName,n_ch_dat,n_bytes))
    return os.path.getsize(FileName)//n_ch_dat//n_bytes
        
 


##################################
######## progress bar stuff ###############
##################################

class sample_counter_widget:
    samples = 0
    def update(self,pbar):
        return str(self.samples)
    
class spike_counter_widget:
    spikes = 0
    def update(self,pbar):
        return str(self.spikes)
        
def fixed_samples_prog(n_samples):
    spike_counter = spike_counter_widget()
    pbar = ProgressBar(widgets=[SimpleProgress()," (",Percentage()," )", " samples processed. ",Bar(">")," ",spike_counter, " spikes found."],maxval=n_samples).start()
    return pbar,spike_counter
    
def fixed_spikes_prog(n_spikes):
    sample_counter = sample_counter_widget()
    pbar = ProgressBar(widgets=[sample_counter, " samples processed. ",Bar(">"),SimpleProgress()," (",Percentage()," )", " spikes found"],maxval=n_spikes).start()
    return pbar,sample_counter

class spikes_and_samples_bar:

    def __init__(self,total_spikes,total_samples):
        if total_spikes is not None:
            self.pbar,self.other_counter = fixed_spikes_prog(total_spikes)
            self.update = self.fixed_spikes_update
        elif total_samples is not None:
            self.pbar,self.other_counter = fixed_samples_prog(total_samples)
            self.update = self.fixed_samples_update

    def fixed_samples_update(self,spikes,samples):
        self.other_counter.spikes = spikes
        self.pbar.update(samples)
        
    def fixed_spikes_update(self,spikes,samples):
        self.other_counter.samples = samples
        self.pbar.update(spikes)
        
    def finish(self):
        self.pbar.finish()


                        
###################################
######## Test Functions ###########
###################################


#################
### for pbar: ###
#################

def test1():
    import time
    
    fixed_spikes,sc = fixed_spikes_prog(100)
    for i in range(100):
        sc.samples = 2*i
        fixed_spikes.update(i)
        time.sleep(.03)    

    fixed_samples,sc = fixed_samples_prog(100)
    for i in range(100):
        sc.spikes = 2*i
        fixed_samples.update(i)
        time.sleep(.03)
        
def test2():
    import time
    
    bar = spikes_and_samples_bar(100,None)
    for i in range(100):
        bar.update(i,i)
        time.sleep(.1)
    bar.finish()
    
    bar = spikes_and_samples_bar(None,100)
    for i in range(100):
        bar.update(i,i)
        time.sleep(.1)
    bar.finish()
