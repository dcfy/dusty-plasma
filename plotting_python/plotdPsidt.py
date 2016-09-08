import gk1
import matplotlib.pyplot as plt
import numpy as np

def smooth_1d(f, width=11, window='hanning', length='origin'):
    '''
    Data smoothing functions adpated from http://www.scipy.org/Cookbook/SignalSmooth
    "window" can be one of flat, bartlett, blackman, hamming, hanning, kaiser
    '''
    if width<3:
        return f
    s=np.r_[f[width-1:0:-1],f,f[-1:-width:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(width)')
    _f_=np.convolve(w/w.sum(),s,mode='valid')

    if length == 'origin':
        #return _f_[width/2-1:-width/2]
        r = len(_f_) - len(f) # r is either odd or even
        return _f_[r/2:-r/2]
    elif length == 'full':
        return _f_

"""
fname = 'gem_ByFlux_20.h5'
byFlux_file = gk1.DynVectorFile(fname)
time = byFlux_file.getTime()
byFlux = byFlux_file.getData()
byFlux_file.close()

plt.plot(time, byFlux, color='k')
plt.show()
"""

fname_fmt = 'gem_ByFlux_%d.h5'
fnames = []
for frame in range(1,20):
    fnames.append(fname_fmt%frame)
time, data = gk1.merge_dyn_vector(fnames)
# larger width means smoother data
data = smooth_1d(data, width=len(data)/20)
#print len(data)/20
ddt = 0.5*(data[1:] - data[:-1])/(time[1:] - data[:-1])

plt.figure()
plt.plot(time, data, color='k', lw=2)
plt.title('Psi')
plt.xlabel('time')
plt.ylabel('ByFlux')

plt.figure()
plt.plot(time[:-1], ddt, lw=2)
tmax = time[np.argmax(ddt)]
ddtmax = np.amax(ddt)
plt.axvline(tmax, linestyle='dashed')
plt.title("dPsi/dt=%g; max at $ t = %g$"%(ddtmax,tmax))
plt.xlabel('time')
plt.ylabel('dByFlux/dt')

plt.show()
