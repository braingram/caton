##########################################################
#### These are the parameters you might want to change ###
##########################################################

T_BEFORE = .0005 # time before peak in extracted spike
T_AFTER = .0005 # time after peak in extracted spike
T_JOIN_CC = .0005 # maximum time between two samples for them to be "contiguous" in detection step
F_LOW = 500. # low pass frequency (Hz)
THRESH_SD = 4.5 # threshold for detection. standard deviations of signal
DETECT_POSITIVE = True # detect spikes with positive threshold crossing
DTYPE = "i2" # ">i2" (> means big-endian), "i4", "f2"
# see http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html#arrays-dtypes-constructing
SEPARATE_CHANNELS_PCA = True # Extract features on channels separately. 
#If False, then I lump all channels together and do PCA, and use that for feature extraction
REGET_FEATURES = False # Extract features from current spikes using PCA. Otherwise use pre-computed waveforms.
SORT_CLUS_BY_CHANNEL = False # Sort clusters by the channel where the peak occurs
# Otherwise clusters are sorted high amplitude to low amplitude
FPC = 3 # Features per channel
BUTTER_ORDER = 3 # Order of butterworth filter
CHUNK_SIZE = 20000   # number of time samples used in chunk for filtering and detection
CHUNK_OVERLAP = 200 # number of samples that chunks overlap in time
CHUNKS_FOR_THRESH = 5 # number of chunks used to determine threshold for detection
INTERP_METHOD = 'linear' # "linear' or 'cubic'. interpolation method when extracting spikes after detection (since peak isn't exactly at sample). Linear is faster. Note that --fast uses no interpolation at all.
