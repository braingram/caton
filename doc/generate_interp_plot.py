import sys,os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath('__file__'))))

from caton.interp_stuff import *
import numpy as np
import matplotlib.pyplot as plt

low = 0
hi = 20
smallstep = .05
bigstep = 1
peak = 9.4
f = lambda x: -np.exp(-(x-peak)**2)

x =  np.arange(low,hi,smallstep)
y = f(x)
x -= x[y.argmin()]

x1 = np.arange(low,hi,bigstep)
y1 = f(x1)
x1 -= x1[y1.argmin()]


x_interp = {}
y_interp = {}
for kind in ['linear','quadratic','cubic']:
    x_interp[kind] = x1.copy()
    y_interp[kind] = interp_around_peak(y1.copy().reshape(-1,1),y1.argmin(),0,10,10,pad=True,kind=kind).flatten()
    x_interp[kind] -= x_interp[kind][y_interp[kind].argmin()]


fig = plt.figure(figsize=(9,9))

ax1 = fig.add_subplot(2,1,1)
ax1.set_title("Waveforms after alignment")
ax1.plot(x,y,'black')
ax1.plot(x1,y1,'blue')

ax2 = fig.add_subplot(2,1,2,sharex=ax1)
ax2.set_title("Errors")
ax2.axhline(0,color='black')
ax2.vlines(x1,[0],y1-y[np.searchsorted(x,x1)],color='blue',lw=2)

for kind,color,offset in zip(['linear','quadratic','cubic'],['red','orange','green'],[.15,.3,.45]):
    ax1.plot(x_interp[kind],y_interp[kind],color=color)
    ax2.vlines(x_interp[kind]+offset,[0],y_interp[kind]-y[np.searchsorted(x,x_interp[kind])],color=color,lw=4)
ax1.legend(("original function","no interpolation","linear splines","quadratic splines","cubic splines"))

plt.savefig("generated_data/interpolation.png")