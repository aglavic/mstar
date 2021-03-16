#-*- coding: utf-8 -*-
'''
Some help functions to analyze and export the data collected in the McStas simulation.
'''

from numpy import *
from copy import deepcopy

def calc_brilliance_transfer(data, data_scale, reference, reference_scale):
  lmin, lmax=map(float, data.info['xlimits'].split())
  l=linspace(lmin, lmax, data.data.shape[0])
  dy=data.data*data_scale
  de=data.errors*data_scale
  ry=reference.data*reference_scale
  re=reference.errors*reference_scale

  bt_y=dy/ry
  bt_e=sqrt((de/ry)**2+(re*dy/ry**2)**2)
  return l, bt_y, bt_e

# Detector parameters for efficiency calculation
DET_INC=5.0*pi/180. #
DET_R1 = 3.08 #µm
DET_R2 = 1.38 #µm
DET_R1a=3.97 #µm
DET_R2a=1.61 #µm
DET_Fa=0.06

def sigma(wavelength):
  return ((2.3*6.022e23/10.33)*0.8*(3844.*wavelength/1.8)*1e-28)/sin(DET_INC) 

def eff(wavelength, R1, R2):
  return (1.-1./(2.*R1*sigma(wavelength))-1./(2.*R2*sigma(wavelength)))+\
         ((exp(-sigma(wavelength)*R2))/(2.*sigma(wavelength)*R2))+\
         ((exp(-sigma(wavelength)*R1))/(2.*sigma(wavelength)*R1))

def det_eff(wavelength):
  return (1.-DET_Fa)*eff(wavelength, DET_R1, DET_R2)+DET_Fa*eff(wavelength, DET_R1a, DET_R2a)

def apply_efficiency(ds):
  dnew=deepcopy(ds)
  I0=dnew.data['p']
  l=dnew.data['L']
  Inew=det_eff(l)*I0
  dnew.data['p']=Inew
  return dnew


################ Reflectivity calculation stuff #################
h=6.6260693E-31 # plancks constant [g·m^2/s]
m_n=1.67493E-24 # Neutron mass [g]
vn_ln=h/m_n*1e10 # conversion factor for velocity in m/s to neutron wavelength in Å
beamline_length=21. # [m] source-detector distance
detector_distance=3.0
lambda_min=2.0
source_frequency=15.0
t_min=lambda_min/(vn_ln/beamline_length)#-(0.0015-0.00175)
t0=0.0002+t_min # 0.00175+t_min #0.00175
l_code='((t-%f)%%%%(%%f/%f)+%f)*%f'%(t0, source_frequency, t_min, vn_ln/beamline_length)
l_code2='((t-%f)+%f)*%f'%(t0, t_min, vn_ln/beamline_length) # for comparing time with wavelength to separate pulses
# The conversion for 2 single and 1 double bandwith pulses to wavelength
l_code_21='where(t<(%(t0)f+2*%(tframe)f), (t-%(t0)f)%%%(tframe)f+%(t_min)f, t-2*%(tframe)f)*%(conversion)f'%{
              't0':t0, 'tframe':1.0/source_frequency, 'conversion': vn_ln/beamline_length, 't_min':t_min}
q_code='4.*pi/lamda*sin(%%f*pi/180.+arctan2(-x,%f))'%(detector_distance)
norm_code='sin(%%f*pi/180.+arctan2(-x,%f))/sin(%%f*pi/180.+arctan2(-x,%f))'%(
            detector_distance, detector_distance)

def calcR(ds, dr, omega, omega_ref, qres=0.01, qmin=0.005, qmax=0.35, mindq=1e-4,
          use_tof=True, skip_pulses=0, lambda_min=2.4, detcorr=False,
          t_overlap=0.001, crop_overlap=True, complex_skip=False):
  res_base=1.0+qres
  q_bins=[qmin]
  while q_bins[-1]<qmax:
    nextq=max(q_bins[-1]*res_base, q_bins[-1]+mindq)
    q_bins.append(nextq)
  q=array(q_bins)


  if complex_skip:
    lc=l_code_21
  elif use_tof:
    lc=l_code%(skip_pulses+1)
  else:
    lc='L'
  if crop_overlap:
    ofltr='(lamda<=%f)&'%eval(l_code%(skip_pulses+1), {'t': (t0+(skip_pulses+1)/source_frequency-t_overlap)})
  else:
    ofltr=''

  if detcorr:
    ds=apply_efficiency(ds)
    dr=apply_efficiency(dr)

  # generate data binned in q
  ignore, Is=ds.project1d('q', newcols=[('lamda', lc), ('q', q_code%omega)],
                          bins=q, fltr=ofltr+'(lamda>=%f)'%lambda_min)
  ignore, Ir=dr.project1d('q', newcols=[('lamda', lc), ('q', q_code%omega)],
                          bins=q, fltr=ofltr+'(lamda>=%f)'%lambda_min,
                          norm=norm_code%(omega, omega_ref))
  return (q[:-1]+q[1:])/2., Is, Is/Ir

def calcStat(ds, window=None, detcorr=True):
  '''
  Calculate event statistics from a dataset.
  The output is a dictionary of event counts:
    'total':    avg. total cps
    'peak':     peak total cps
    'roi':      avg. cps per mm² for ROI on detector
    'roi_peak': peak cps per mm² for ROI on detector
    'l':        array of wavelength values
    'Iroi':     array of ROI instantaneous cps per mm² values
    'Itot':     array of total instantaneous cps values
  '''
  if detcorr:
    ds=apply_efficiency(ds)

  if window is None:
    # find detector region with highest count rate
    x, y, I=ds.project2d('x', 'y', bins=20, fltr='x>=-0.06')
    X, Y=meshgrid((x[:-1]+x[1:])/2., (y[:-1]+y[1:])/2.)
    x0=X.flatten()[I.argmax()]
    y0=Y.flatten()[I.argmax()]
    window=[x0-0.005, x0+0.005, y0-0.005, y0+0.005]

  out={}
  l, I=ds.project1d('L', bins=50, fltr='x>=-0.06')
  li=(l[:-1]+l[1:])/2.

  out['total']=I.sum()
  I*=50.
  out['peak']=I.max()
  out['l_peak']=li[I.argmax()]
  out['l']=li
  out['Itot']=I

  l, I=ds.project1d('L', bins=l,
                    fltr='(x>=%g)&(x<=%g)&(y>=%g)&(y<=%g)'%tuple(window))
  I/=1e6*(window[1]-window[0])*(window[3]-window[2])
  out['roi']=I.sum()
  I*=50.
  out['roi_peak']=I.max()
  out['Iroi']=I

  return out


def save_w_header(fname, col_list, info, names, units):
  data=array(col_list).T
  fh=open(fname, 'w')
  fh.write('# '+'\n# '.join(info.splitlines()))
  fh.write('\n')
  fh.write('# '+' '.join(['%-24s'%col for col in names]))
  fh.write('\n')
  fh.write('# '+' '.join(['%-24s'%('['+col+']') for col in units]))
  fh.write('\n')

  savetxt(fh, data)

