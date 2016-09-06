import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import gk1

lightSpeed = 1.
mu0 = 1.
n0 = 1.
mi = 1.
qi=1.
me = 1/25.
qe=-1.
wpe0_wce0 = 3.
epsilon0 = 1/lightSpeed*lightSpeed/mu0
B0 = lightSpeed * np.sqrt(mu0 * n0 * me) / wpe0_wce0
v_A = B0/np.sqrt(mu0*n0*mi)
wci0=qi*B0/mi
wpi0 = np.sqrt(n0 * qi*qi / epsilon0 / mi)
di0 = lightSpeed/wpi0
de0 = di0*np.sqrt(me/mi)

Enorm = B0*v_A
print Enorm
iEnorm = 1./Enorm
tnorm = 1./wci0
xnorm = di0

#def plotOhmY(ts, ix=None, species=['e','p','o']):
#  f = gk1.File('gem_q_%d.h5'%ts,
#          num_fluids=3,fluid_names=['e','p','o'],num_moments=5,
#def plotOhmY(ts, ix=None, species=['e','p']):
def plotOhmY(ts, ix=None, iy=None, species=['e','p']):
  f = gk1.File('gem-1010-ds_q_%d.h5'%ts,
          num_fluids=2,fluid_names=['e','p'],num_moments=10,
          me=1./25., qe=-1., mp=1., qp=1., mo=16., qo=1.)
  nx,ny,ncomp = f['StructGridField'].shape
  xlo,ylo = f.LowerBounds()
  xup,yup = f.UpperBounds()
  dx = (xup-xlo)/nx
  dy = (yup-ylo)/ny
  idx = 1./dx
  idy = 1./dy
  if ix == None:
    ix = nx/2
    #print 'ix=', ix
  if iy == None:
    iy = ny/2
    #print 'iy=', iy
  x = f.getCoordinates('x')/xnorm
  #print 'x=', x, len(x)
  y = f.getCoordinates('y')/xnorm
  #print 'y=', y, len(y)
  x_cut = x[ix]
  #print 'x_cut=', x_cut
  y_cut = y[iy]
  #print 'y_cut=', y_cut
  fig, ax = plt.subplots()
  #
  Ez = f.getField('Ez', iy=iy)*iEnorm
  #print 'Ez=', Ez, len(Ez)
  ax.plot(x, Ez, label='$ E_z $',color='k',lw=2,alpha=0.8)
  #
  Bx = f.getField('Bx', iy=iy)
  By = f.getField('By', iy=iy)
  #
  for s in species:
    # -VxB term, i.e., convection electric field
    vx_s = f.getField('vx_'+s, iy=iy)*iEnorm
    vy_s = f.getField('vy_'+s, iy=iy)*iEnorm
    _v_sxB = -vx_s*By + vy_s*Bx
    ax.plot(x, _v_sxB, label='$ -(\\mathbf{v}_%s \\times \\mathbf{B} )|_z$'%(s),lw=2,alpha=0.8)
    # TODO: implement more terms
  title = '$t=%g, x=%g$'%(f.Time()/tnorm,y_cut)
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title = title, fontsize='x-large')
  lgd.get_title().set_fontsize('x-large')
  ax.set_xlabel('$x/d_{i0}$', fontsize='x-large')
  ax.set_ylabel('$E/B_0v_{A0}$', fontsize='x-large')
  ax.set_xlim([-20, 20])
  ax.set_ylim([-0.4, 0.1])
  fig.savefig('ohm_Ez_ycut_%02d_x%g.png'%(ts,y_cut), bbox_inches='tight')
  #plt.show()
  plt.close(fig)
  f.close()

#for frame in [10]:
for frame in range(80):
#for frame in range(57,58):
  print("frame %d"%frame)
  plotOhmY(frame)
