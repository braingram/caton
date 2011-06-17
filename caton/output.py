import numpy as np
from xml.etree.ElementTree import ElementTree,Element,SubElement
from utils_misc import switch_ext
import os.path


def write_clu(clus,filepath):
    """writes cluster cluster assignments to text file readable by klusters and neuroscope.
    input: clus is a 1D or 2D numpy array of integers
    output:
        top line: number of clusters (max cluster)
        next lines: one integer per line"""
    clu_file = open( filepath,'w')
    #header line: number of clusters
    n_clu = clus.max()+1
    clu_file.write( '%i\n'%n_clu)
    #one cluster per line
    np.savetxt(clu_file,np.int16(clus),fmt="%i")
    clu_file.close()
    
    
def read_clu(filepath):
    """skip first line, read the rest into an array and return it"""
    return np.loadtxt(filepath,dtype=np.int32,skiprows=1)

def write_fet(feats,filepath,samples=None):
    """writes array of feature vectors to text file readable by klusters and klustakwik
    FOR KLUSTERS, YOU MUST GIVE samples! OTHERWISE IT WILL CRASH WITHOUT SENSIBLE MESSAGE
    input: feats is a 2D ndarray of floats. n_vectors x n_features per vector
        optionally also input samples (times) vector
    output:
        top line: number of features
        next line: one feature vector per line, as integers. last column is time vector if specified.
        last line ends in newline."""
    feat_file = open(filepath,'w')
    #rescaling features so they line between -16383 and 16384
    #feat_scaling_factor = 16000./max(feats.max(),-feats.min())
    #feats *= feat_scaling_factor
    feats = np.int32(feats)
    if samples is not None:
        feats = np.hstack( (feats,samples.reshape(-1,1)) )
    #header line: number of features
    feat_file.write( '%i\n'%feats.shape[1] )
    #next lines: one feature vector per line
    np.savetxt(feat_file,feats,fmt="%i")
    feat_file.close()
    
def read_fet(filepath):
    """reads feature file and returns it as an array. note that the last
    column might contain the times"""
    #skip first line and read the rest
    return np.loadtxt(filepath,dtype=np.int32,skiprows=1).astype(np.float32)
    
    
def write_res(samples,filepath):
    """input: 1D vector of times shape = (n_times,) or (n_times, 1)
    output: writes .res file, which has integer sample numbers"""
    np.savetxt(filepath,samples,fmt="%i")
    
def read_res(filepath):
    """reads .res file, which is just a list of integer sample numbers"""
    return np.loadtxt( filepath,dtype=np.int32)

def write_spk(waves,filepath,nonzero=None):
    """input: waves: 3D array of waveforms. n_spikes x n_channels x n_samples
    nonzero [optional]: 2D boolean array n_spikes x n_channels
    rescaled to signed 16-bit integer and written to file filedir/filebase.spk.1"""
    #wave_scaling_factor = 16000./max(waves.max(),-waves.min())
    if nonzero is not None:
        waves = waves*nonzero.reshape( nonzero.shape + (1,) )
    #waves *= wave_scaling_factor
    waves = np.int16(waves)
    waves.tofile(filepath)
    
def read_spk(filepath,n_ch,n_s):
    return np.fromfile(filepath,dtype=np.int16).reshape(-1,n_s,n_ch)
    
def write_xml(n_ch,n_samp,n_feat,sample_rate,filepath):
    """makes an xml parameters file so we can look at the data in klusters"""
    parameters = Element('parameters')
    acquisitionSystem = SubElement(parameters,'acquisitionSystem')
    SubElement(acquisitionSystem,'nBits').text = '16'
    SubElement(acquisitionSystem,'nChannels').text = str(n_ch)
    SubElement(acquisitionSystem,'samplingRate').text = str(int(sample_rate))
    SubElement(acquisitionSystem,'voltageRange').text = '20'
    SubElement(acquisitionSystem,'amplification').text = "1000"
    SubElement(acquisitionSystem,'offset').text = "2048"

    channels = SubElement(SubElement(SubElement(parameters,'channelGroups'),'group'),'channels')
    for i_ch in range(n_ch):
        SubElement(channels,'channel').text=str(i_ch)
    
    group = SubElement(SubElement(SubElement(parameters,'spikeDetection'),'channelGroups'),'group')
    channels = SubElement(group,'channels')
    for i_ch in range(n_ch):
        SubElement(channels,'channel').text=str(i_ch)
    SubElement(group,'nSamples').text = str(n_samp)
    SubElement(group,'peakSampleIndex').text = str(n_samp//2)
    SubElement(group,'nFeatures').text = str(n_feat)
    
    indent_xml(parameters)
    ElementTree(parameters).write(filepath)    
    

def indent_xml(elem, level=0):
    """input: elem = root element
    changes text of nodes so resulting xml file is nicely formatted.
    copied from http://effbot.org/zone/element-lib.htm#prettyprint"""
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
            
            
def get_pars_from_xml(xmlpath):
    assert os.path.exists(xmlpath)
    root = ElementTree().parse(xmlpath)
    acquisitionSystem = root.find('acquisitionSystem')
    n_channels = int(acquisitionSystem.find('nChannels').text)
    sample_rate = np.float32(acquisitionSystem.find('samplingRate').text)
    return n_channels,sample_rate

def get_pars_from_xml2(xmlpath):
    root = ElementTree().parse(xmlpath)
    n_channels = int(search_etree(root,"nChannels"))
    sample_rate = float(search_etree(root,"samplingRate"))
    s_total = int(search_etree(root,"nSamples"))
    s_before = int(search_etree(root,"peakSampleIndex"))
    s_after = s_total-s_before
    return n_channels,sample_rate,s_before,s_after

def walk_etree(root):
    yield root.tag,root.text
    for child in root.getchildren():
        for tag,text in walk_etree(child):
            yield tag,text
    
def search_etree(root,the_tag):
    for tag,text in walk_etree(root):
        if tag==the_tag: return text

def get_dat_pars(DatFileName):
    xmlpath = switch_ext(DatFileName,'xml')
    if os.path.exists(xmlpath):
        n_ch_dat,sample_rate = get_pars_from_xml(xmlpath)
    else:
        n_ch_dat,sample_rate = get_pars_from_prompt()
        write_xml(n_ch_dat,0,0,sample_rate,xmlpath)
        print("writing parameters in xml file %s"%xmlpath)
    return n_ch_dat,sample_rate
                                          
def get_pars_from_prompt():
    print("Could not find xml file with parameters.")
    n_ch_dat = input("How many channels in .dat file?\t")
    sample_rate = input("What is the sample rate?\t")
    return n_ch_dat,sample_rate
