#!/usr/bin/env python
#-*- coding: utf-8 -*-
'''
Main script that runs the full analysis for the M-STAR McStas simulation results.
'''

import sys
import os
from numpy import *
from matplotlib import pyplot

from mcstas_reader import McSim
from reduce_help import calcR, l_code, vn_ln, beamline_length

RESULT_PATH="../results/"

def do_focusing():
    resopt=dict(lambda_min=2.4, qres=0.01, mindq=5e-5, qmax=0.35)

    r10=McSim(RESULT_PATH+'reference_10/')['detectorEvents']
    r34=McSim(RESULT_PATH+'reference_34/')['detectorEvents']
    s10=McSim(RESULT_PATH+'sample_10/')['detectorEvents']
    s34=McSim(RESULT_PATH+'sample_34/')['detectorEvents']

    q, I, R=calcR(s10, r10, 1.0, 1.0, **resopt)
    q2, I2, R2=calcR(s34, r34, 3.4, 3.4, **resopt)

    r20=McSim(RESULT_PATH+'reference_20/')['detectorEvents']
    s20=McSim(RESULT_PATH+'sample_20/')['detectorEvents']
    qs, Is, Rs=calcR(s20, r20, 2.0, 2.0, complex_skip=True, **resopt)

    savetxt(RESULT_PATH+'ni_10x10.dat',
            array([q, R, I*15/4, R2, I2*15/4, Rs, Is*15/4]).T,
            fmt="%16e",header='M-STAR focusing reflectometry on 10x10 nickel on silicon\n'+
            ''.join(['%17s'%si for si in
                     ['q','R(1.0)','I(1.0)', 'R(3.4)','I(3.4)','Rskip(2.0)','Iskip(2.0)']])[3:])


    s10=McSim(RESULT_PATH+'sample_2x2_10/')['detectorEvents']
    s34=McSim(RESULT_PATH+'sample_2x2_34/')['detectorEvents']

    q, I, R=calcR(s10, r10, 1.0, 1.0, **resopt)
    q2, I2, R2=calcR(s34, r34, 3.4, 3.4, **resopt)

    savetxt(RESULT_PATH+'ni_2x2.dat',
            array([q, R, I*15/4, R2, I2*15/4]).T,
            fmt="%16e",header='M-STAR focusing reflectometry on 2x2 nickel on silicon\n'+
            ''.join(['%17s'%si for si in ['q','R(1.0)','I(1.0)', 'R(3.4)','I(3.4)']])[3:])

def do_collimated(lmin=3.9, lmax=16.0):
    r=McSim(RESULT_PATH+'reference_coll/')['detectorEvents']

    ledges=1.01**(arange(log(lmin)/log(1.01), log(lmax)/log(1.01)+2, 1)-0.5)
    l=1.01**(arange(log(lmin)/log(1.01), log(lmax)/log(1.01)+1, 1))

    source_frequency=15.0
    t_min=lmin/(vn_ln/beamline_length)  # -(0.0015-0.00175)
    t0=0.0002+t_min  # 0.00175+t_min #0.00175
    l_code='((t-%f)%%%%(%%f/%f)+%f)*%f'%(t0, source_frequency, t_min, vn_ln/beamline_length)
    ignore, ref=r.project1d('lamda', bins=ledges, newcols=[('lamda', l_code%1)])

    for omega in [0.4, 1.5, 6.0]:
        s=McSim(RESULT_PATH+'sample_coll_%02i/'%int(omega*10))['detectorEvents']
        BS=3.0*tan(omega*pi/180*0.02)
        ignore, I=s.project1d('lamda', bins=ledges, newcols=[('lamda', l_code%1)],
                              fltr='abs(x)<%f'%(BS*3))
        q=4.*pi/l*sin(omega*pi/180.)

        savetxt(RESULT_PATH+'ni_coll_%02i.dat'%int(omega*10),
                array([q, I/ref/omega**2, I*15/4]).T,
                fmt="%16e",header='M-STAR collimated reflectometry on 10x10 nickel on silicon\n'+
                ''.join(['%17s'%si for si in ['q','R(1.0)','I(1.0)']])[3:])


if __name__=="__main__":
    do_focusing()
    do_collimated()
