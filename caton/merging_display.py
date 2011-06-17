from __future__ import division
import numpy as np, matplotlib.pyplot as plt, itertools as it

NO_VALUE = -423

MIN_COUNT_ARROW = 10
MIN_FRAC_ARROW = .03
MAX_ROWS = 40

def row_dtype(n_ch):
    return np.dtype([("rank",int),("left",int),("good",bool),("ch_amp",np.float32,n_ch)])

def arrow(x1,y1,x2,y2,width):
    plt.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(edgecolor="blue",width=width,headwidth=2*width))

    
def merge_diagnostics(n_ch,key2subset,key2rank,key2left,key2good,key2spkmean,fromto2stolen):
    n_clu = len(key2subset)
    keys = key2subset.keys()
    clu_data = np.recarray(n_clu,dtype=row_dtype(n_ch))
    key2rowind = {}
    for i_row,key in enumerate(sorted(keys,key=lambda key: np.mean(key2subset[key]))):
        key2rowind[key] = i_row
        ch_amp = np.empty(n_ch,dtype=np.float32)
        ch_amp.fill(NO_VALUE)
        ch_amp[key2subset[key]] = key2spkmean[key].ptp(axis=0)
        clu_data[i_row] = (key2rank[key],key2left[key],key2good[key],ch_amp)
        
        
    stolen_arr = np.zeros((n_clu,n_clu),dtype=np.int32)
    for ((fromkey,tokey),count) in fromto2stolen.items(): stolen_arr[key2rowind[fromkey],key2rowind[tokey]] = count    
    plot_clu_data(clu_data[:MAX_ROWS],stolen_arr[:MAX_ROWS,:MAX_ROWS])
    
def text(x,y,v,**kw):
    textkw = dict(horizontalalignment="center",verticalalignment="bottom")
    textkw.update(kw)
    if v==NO_VALUE: s = ""
    elif str(v).find(".") != -1: s = "%.2f"%v
    else: s = str(v)
    plt.text(x,y,s,**textkw)
    
def plot_clu_data(clu_data,stolen_arr):    
    n_ch = clu_data.ch_amp.shape[1]
    y_top = .9
    x_rank,x_left = .05,.1
    x_ch_amp = np.linspace(.15,.9,n_ch,endpoint=False)
    y_rows = np.linspace(.9,.1,clu_data.shape[0],endpoint=False)
    dx,dy = x_ch_amp[1]-x_ch_amp[0],y_rows[1]-y_rows[0]
    for i_ch,x_ch in enumerate(x_ch_amp): text(x_ch+dx/2,y_top,i_ch)
    ax = plt.gca()
    text(x_rank,y_top,"rank")
    text(x_left,y_top,"spk left")
    
    rescale = lambda amp: 256/np.log(clu_data.ch_amp.max())*np.log(amp)    
    def my_cmap(x,b):
        if x == NO_VALUE:
            return (.7,1,.7) if b else (1,.7,.7)
        else: return plt.cm.gist_yarg(int(rescale(x)))
    
    for (row,y_row) in zip(clu_data,y_rows):        
        text(x_rank,y_row+dy,row.rank)
        text(x_left,y_row+dy,row.left)
        for x_ch,ch_amp in zip(x_ch_amp,row.ch_amp): 
            r = plt.Rectangle((x_ch,y_row),dx,dy,fc=my_cmap(ch_amp,row.good))
            ax.add_patch(r)    
    x_line = 0
    ongrid = lambda x: x*(x_ch_amp[-1]-x_ch_amp[0])+x_ch_amp[0]
    for (i_from,i_to), count in np.ndenumerate(stolen_arr):
        stole_frac = count/clu_data.left[i_from]
        if count >= MIN_COUNT_ARROW and stole_frac > MIN_FRAC_ARROW:
            arrow(ongrid(x_line),y_rows[i_from]+dy/2,ongrid(x_line),y_rows[i_to]+dy/2,7*stole_frac)
            x_line = (x_line + .03)%1
    plt.show()
    
def isize(set1,set2):
    return len(set1.intersection(set2))
    
if __name__ == "__main__":
    data = np.recarray(3,dtype=row_dtype(2))
    plot_clu_data(data)