#!/usr/bin/env python
from __future__ import division
import numpy as np
from utils_graphs import contig_segs
from utils_misc import *
import tableprint
import output
import os


    
def print_isect_table_from_dirs(sad_dir,happy_dir):
    rc_path = find_file_with_ext(sad_dir,'.clu.1',ex_if_not_found = True)
    rt_path = find_file_with_ext(sad_dir,'.res.1',ex_if_not_found = True)
    mc_path = find_file_with_ext(happy_dir,'.clu.1',ex_if_not_found = True)
    mt_path = find_file_with_ext(happy_dir,'.res.1',ex_if_not_found = True)        

    
    xml_path1 = find_file_with_ext(happy_dir,'.xml')
    xml_path2 = find_file_with_ext(sad_dir,'.xml')
    xml_path = xml_path1 if xml_path1 else xml_path2    
    
    
    oprint("Comparing cluster file %s with time file %s to cluster file %s with time file %s"
          %map(os.path.basename,(rc_path,rt_path,mc_path,mt_path)))
    from xml.etree.ElementTree import ElementTree
    sample_rate = float(ElementTree().parse(xml_path).find('acquisitionSystem').find('samplingRate').text)
    print_isect_table(output.read_res(rt_path),output.read_clu(rc_path),
                       output.read_res(mt_path),output.read_clu(mc_path),sample_rate)

def print_report_from_dirs(extra_dir,intra_dir):
    extra_clu_path = find_file_with_ext(extra_dir,'.clu.1',ex_if_not_found = True)
    extra_res_path = find_file_with_ext(extra_dir,'.res.1',ex_if_not_found = True)
    intra_res_path = find_file_with_ext(intra_dir,'.res.1',ex_if_not_found = True)        

    
    xml_path = find_file_with_ext(extra_dir,'.xml')
    
    
    oprint("Comparing extracellular units from res file %s, clu file %s to intracellular unit with time file %s"
          %tuple(map(os.path.basename,(extra_clu_path,extra_res_path,intra_res_path))))
    from xml.etree.ElementTree import ElementTree
    sample_rate = float(ElementTree().parse(xml_path).find('acquisitionSystem').find('samplingRate').text)
    
    meas_clu_counts,isect_counts,n_undetected = match_real_to_meas(output.read_res(intra_res_path),output.read_res(extra_res_path),
                                                                   output.read_clu(extra_clu_path),sample_rate)
    print_report(meas_clu_counts,isect_counts,n_undetected,len(output.read_res(intra_res_path)))
    
def match_real_to_meas(real_times,meas_times,meas_clus,sample_rate):
    TOLERANCE_TIME = .001
    t_tol_half = int(TOLERANCE_TIME*sample_rate/2)    
    all_times = np.zeros(max(real_times[-1],meas_times[-1])+1000,dtype=np.int16)
    for t in real_times:
        all_times[t-t_tol_half:t+t_tol_half] = 1
    n_meas_clus = meas_clus.max()+1
    isect_counts = np.zeros(n_meas_clus,dtype=np.int32)
    for clu,t in zip(meas_clus,meas_times):
        if all_times[t] != 0:
            isect_counts[clu]+=1
            all_times[t] = -1
    undetected = 0
    for t in real_times:
        if (all_times[t-t_tol_half:t+t_tol_half] == 1).all():
            undetected+=1
    return np.bincount(meas_clus),isect_counts,undetected

def print_report(meas_clu_counts,isect_counts,n_undetected,n_intra):
    best_clu = isect_counts.argmax()
    oprint("Total intra: %i"%n_intra)
    oprint("Undetected: %i"%n_undetected)
    n_good = isect_counts[best_clu]
    n_wrong = meas_clu_counts[best_clu]-n_good
    oprint("Wrong spikes in best cluster: %i ==> Type I error rate: %.1f%%"%(n_wrong,100*n_wrong/n_intra))
    n_missing = n_intra-n_good
    oprint("Intra spikes missing from best cluster: %i ==> Type II error rate: %.1f%%"%(n_missing,100*n_missing/n_intra))
    

            
def print_isect_table(sad_times,sad_clus,happy_times,happy_clus,sample_rate):
    TOLERANCE_TIME = .001
    tolerance_samples = int(TOLERANCE_TIME*sample_rate)
    
    sad_clu_labels, = np.bincount(sad_clus).nonzero()
    happy_clu_labels, = np.bincount(happy_clus).nonzero()
    
    max_sad_label = sad_clus.max()
    max_happy_label = happy_clus.max()
    
    intersection_table = np.zeros((max_happy_label+1,max_sad_label+1),dtype=np.int32)

    closest_happy_to_sad = closest_inds_unsorted(happy_times,sad_times)
    happy_intersect_none = np.bincount(happy_clus)
    sad_intersect_none = np.zeros(max_sad_label+1,dtype=np.int32)
    
    for sad_ind,sad_time in enumerate(sad_times):
        sad_clu = sad_clus[sad_ind]
        happy_clu = happy_clus[closest_happy_to_sad[sad_ind]]
        if np.abs(happy_times[closest_happy_to_sad[sad_ind]]-sad_time)<=tolerance_samples:
            happy_intersect_none[happy_clu] -= 1
            intersection_table[happy_clu,sad_clu] += 1
        else:
            sad_intersect_none[sad_clu] += 1        
            
            
    results_table = np.zeros((len(happy_clu_labels)+3,len(sad_clu_labels)+3),dtype='|S8')
    results_table[1:-2,1:-2] = intersection_table[np.ix_(happy_clu_labels,sad_clu_labels)]
    results_table[1:-2,-2] = happy_intersect_none[happy_clu_labels]
    results_table[-2,1:-2,] = sad_intersect_none[sad_clu_labels]
    results_table[1:-2,0] = happy_clu_labels
    results_table[0,1:-2] = sad_clu_labels
    results_table[0,-2] = "None"
    results_table[-2,0] = "None"
    results_table[1:-2,-1] = np.bincount(happy_clus)[happy_clu_labels]
    results_table[-1,1:-2] = np.bincount(sad_clus)[sad_clu_labels]
    results_table[0,-1] = "Total"
    results_table[-1,0] = "Total"
    results_table[0,0] = '\\'
    
    
    if results_table.shape[1] > results_table.shape[0]:
        results_table = results_table.T
    
    oprint(tableprint.indent(list(results_table), 
                    hasHeader=True, separateRows=True,
                     prefix='| ', postfix=' |'))
    


def closest_inds_sorted(base,inserted):
    i_base_after_ins = np.searchsorted(base,inserted)
    i_base_before_ins = np.searchsorted(base,inserted)-1
    base2 = np.r_[-1000000000,base,1000000000]
    after_closer = base2[i_base_after_ins+1]-inserted < inserted-base2[i_base_before_ins+1]
    return after_closer*(i_base_after_ins)+(1-after_closer)*(i_base_before_ins)

def closest_inds_unsorted(base,inserted):
    "output y_i = index of base closest to inserted_i"
    sbase = np.sort(base)
    sins = np.sort(inserted)
    sclose = closest_inds_sorted(sbase,sins)
    
    baseorder = np.argsort(base)
    insrank = np.argsort(np.argsort(inserted))
    #unsort domain and range
    return baseorder[sclose[insrank]]
    
                
def test_match():
    real_times = np.sort(np.random.randint(2400020,size=1000))
    real_clus = np.random.randint(3,size=real_times.size)
    meas_times = np.sort(np.random.randint(2400020,size=9000))
    meas_clus = np.random.randint(4,size=meas_times.size)
    meas_clus += 2
    sample_rate = 10000.
    print_isect_table(real_times,real_clus,meas_times,meas_clus,sample_rate)

def test_report():
    intra_dir = "/home/joschu/Data/d11221/d11221.002_intracellular"
    extra_dir = "/home/joschu/Data/d11221/d11221.002_tet_batch"
    print_report_from_dirs(extra_dir,intra_dir)
    
def main():
    from optparse import OptionParser
    usage = "%prog intra_directory extra_directory"
    parser = OptionParser(usage)    
    (opts,args) = parser.parse_args()
    intra_dir,extra_dir = args
    print_report_from_dirs(intra_dir,extra_dir)
    
if __name__ == '__main__':
    test_report()
    #main()
