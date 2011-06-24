from __future__ import division, with_statement
import os,numpy as np, tempfile,collections, itertools as it
from output import write_fet,read_clu
from CEM_extensions import subset_inds
from merging_display import merge_diagnostics
from utils_misc import time_fun
import probe_stuff
from os.path import join

DEBUG = False
MINCLUSTERS = 3
MAXCLUSTERS = 14
ACCEPTABLE_FRAC = .8
MIN_CLU_SIZE = 10

    
def spike_subsets(ST_nc,ChSubsets):
    ch2spikes = [np.flatnonzero(ST_n) for ST_n in ST_nc.transpose()]
    return [reduce(np.union1d,
                   [ch2spikes[ch] for ch in subset]) 
            for subset in ChSubsets]
    
def clu_mag(Mean_mf):
    return (Mean_mf**2).sum(axis=1)


def cluster_withsubsets(spike_table,clusterdir,reorder_clus=True):
    "TODO: write docstring"
    
    if reorder_clus: print "Cluster reordering not implemented!"
    ST_nc = np.bool8(spike_table.cols.st[:])
    Fet_nc3 = spike_table.cols.fet[:]    
    
    ChSubsets = probe_stuff.SORT_GROUPS
    SpkSubsets = spike_subsets(ST_nc,ChSubsets)    
    print("%i subsets total"%len(SpkSubsets))
    n_spikes,n_ch,_FPC = Fet_nc3.shape
    
    key2subset, key2members, key2spkmean, key2mag = {},{},{},{}
    for i_subset,ChHere,SpkHere in zip(it.count(),ChSubsets,SpkSubsets):        
        print("Sorting channels %s"%ChHere.__repr__())
        FetHere_nc3 = Fet_nc3[np.ix_(SpkHere,ChHere)] # features of spikes in this subset
        CluArr = klustakwik_cluster(FetHere_nc3,'/'.join((clusterdir,"cluster_%i" % i_subset)))
        CluMembersList = [(SpkHere[inds]) for inds in subset_inds(CluArr)] #go back to original indices
        # We are ignoring cluster 0 here, because of [1:] above. No not now
        for (i_clu,Members) in enumerate(CluMembersList):
            if len(Members) > MIN_CLU_SIZE:
                SpkMean = np.array([spike_table[member]["wave"][:,ChHere] for member in Members]).mean(axis=0)
                key = (i_subset,i_clu)
                key2subset[key]=ChHere
                key2members[key] = Members
                key2spkmean[key] = SpkMean
                key2mag[key] = SpkMean.ptp(axis=0).sum()
        
    ImprovingKeys = sorted(key2mag.keys(),key = lambda key: key2mag[key])    
    #problem: most spikes aren't members of any cluster?!
    
    key2oldcount = dict((key,len(members)) for key,members in key2members.items())
    FinalClu = np.zeros(n_spikes,dtype=np.dtype([("subset",int),("clu",int)]))

    # maybe i should have a key2int kind of function?
    fromto2stolen = collections.defaultdict(int)
    for key in ImprovingKeys:
        if DEBUG: 
            for oldkey in FinalClu[key2members[key]]: fromto2stolen[tuple(oldkey),key] += 1
        FinalClu[key2members[key]] = key
    for fromkey,tokey in fromto2stolen.keys(): 
        if DEBUG:
            if fromkey == (0,0): del fromto2stolen[(fromkey,tokey)]
        
    key2newcount = dict((key,((FinalClu["subset"] == key[0]) & (FinalClu["clu"] == key[1])).sum()) for key in ImprovingKeys)    
    key2good = dict((key,
                     key2newcount[key]/key2oldcount[key] > ACCEPTABLE_FRAC and
                     key2oldcount[key] > MIN_CLU_SIZE)
                    for key in ImprovingKeys)

    good_keys = filter(lambda key: key2good[key],reversed(ImprovingKeys))
    
    #with open("counts.txt","w") as fd:
    #    for i_clu,(new,old) in enumerate(zip(NewCount,OrigCount)):
    #        fd.write("%i: %i/%i\n"%(i_clu,new,old) if new/old < .8 else "%i: %i/%i ==> %i\n"%(i_clu,new,old,RelabelArr[i_clu]))

    # problem: relabel cluster indices so they're in the right order
    
    key2rank = dict((key,rank) for (rank,key) in enumerate(reversed(ImprovingKeys)))
    key2left = dict((key,len(members)) for key,members in key2members.items())
    
    if DEBUG: 
        merge_diagnostics(n_ch,key2subset,key2rank,key2left,key2good,key2spkmean,fromto2stolen)
    key2ind = dict((key,ind) for (ind,key) in enumerate(sorted(good_keys,key=lambda key: np.mean(key2subset[key]))))
    FinalCluInd = np.array([key2ind.get(tuple(key),0) for key in FinalClu],dtype=np.int32)
    return FinalCluInd


#def reorder_clus(Fet_nc3,Clu_n,fixed = [0]):
    #print("Reordering clusters")
    #n_spikes,n_ch,fpc = Fet_nc3.shape
    #n_clu = Clu_n.max()+1
    #Mean_mf = class_means(to2d(Fet_nc3).astype(np.float32),Clu_n,n_clu)   
    #Mean_mc3 = Mean_mf.reshape(n_clu,n_ch,fpc)
    #PeakCh_m = (Mean_mc3**2).sum(axis=2).argmax(axis=1)
    
    #ClusToReorder = np.arange(n_clu,dtype=np.int32)
    #ClusToReorder = np.delete(ClusToReorder,fixed)
    #Old2New = np.empty(n_clu,dtype=np.int32)
    #Old2New[fixed] = fixed
    #Old2New[ClusToReorder] = ClusToReorder[np.argsort(PeakCh_m[ClusToReorder])]
    
    #return Old2New[Clu_n]

@time_fun
def klustakwik_cluster(Fet_nc3, clusterdir):
    kk_path = "KlustaKwik"
    if not (os.path.exists(clusterdir)):
        os.makedirs(clusterdir)
    kk_input_filepath = join(clusterdir,'k_input.fet.1')
    kk_output_filepath = join(clusterdir,'k_input.clu.1')
    Fet_nf = Fet_nc3.reshape(len(Fet_nc3),-1)
    write_fet(Fet_nf,kk_input_filepath)
    n_fet = Fet_nf.shape[1]
    
    os.system( ' '.join([kk_path, kk_input_filepath[:-6], '1', '-UseFeatures', '1'*n_fet,'-MinClusters',str(MINCLUSTERS),'-MaxClusters',str(MAXCLUSTERS),'-Screen','0']) )
    return read_clu(kk_output_filepath)

    
