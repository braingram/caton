import numpy as np
from scipy.signal import cspline1d,cspline1d_eval
from scipy.interpolate import interp1d

def abc(x_3,y_3):
    M = np.vstack((x_3**2,x_3,np.ones_like(x_3)))
    return np.linalg.solve(M.T,y_3)

def max_t(a_b_c):
    return -a_b_c[1]/(2*a_b_c[0])

def interp_around(X_sc,s_fracpeak,s_before,s_after,kind='cubic'):
    n_c = X_sc.shape[1]
    n_s = s_before+s_after
    Out_sc = np.empty((n_s,n_c),dtype=np.float32)
    for i_c in xrange(n_c):
        if kind == 'cubic':
            coeffs = cspline1d(X_sc[:,i_c])
            Out_sc[:,i_c] = cspline1d_eval(coeffs,
                                       newx=np.arange(s_fracpeak - s_before,s_fracpeak+s_after,dtype=np.float32))
        elif kind == "linear":
            Out_sc[:,i_c] = interp1d(np.arange(X_sc.shape[0]),X_sc[:,i_c],
                                     bounds_error=True,kind=kind)(np.arange(s_fracpeak - s_before,s_fracpeak+s_after,dtype=np.float32))
        else: raise Exception("kind must be 'linear' or 'cubic'")
    return Out_sc

def interp_around_peak(X_sc,i_intpeak,c_peak,s_before,s_after,pad=False,kind='cubic'):    

    a_b_c = abc(np.arange(i_intpeak-1,i_intpeak+2,dtype=np.float32),
                X_sc[i_intpeak-1:i_intpeak+2,c_peak])
    s_fracpeak = max_t(a_b_c)
    if pad:
        return interp_around(np.vstack((X_sc[0],X_sc,X_sc[-1])),s_fracpeak+1,s_before,s_after,kind=kind)
    else:
        return interp_around(X_sc,s_fracpeak,s_before,s_after,kind=kind)      