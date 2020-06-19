'''Originally written by Mahdi Bahrampouri at Virginia Tech'''
'''Modified by Som Dhulipala (06/18/20)'''


import os
os.chdir('/Users/som/Dropbox/SiteResp_GMSelect/Python_SourceCode')
import pdb
import csv
import EquivLin
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
import degrcurv as dc
error = 0.01

ss = os.getcwd()
def response_fourier(acc_g,dt,w,zeta):
    f = np.fft.fftfreq(d = dt,n = 2**14)
    FA = np.fft.fft(acc_g,n=2**14)
    zeta=0.05
    w = np.float64(np.ravel(w))
    z = np.float64(zeta)
    ii=0
    SA2 = np.ones(len(w))
    for w1 in w:
        fn = w1/(np.pi*2)
        hf = -fn**2/(f**2-fn**2-2*1j*z*f*fn)
        acc_fr_s = FA*hf
        acc_t_s = np.fft.irfft(acc_fr_s)
        SA2[ii] = np.max(np.absolute(acc_t_s))
        ii=ii+1
    return SA2,w

# T = np.array([8.2325])
T=np.array([0.01, 0.025, 0.05,0.075, 0.1, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6,0.7,0.85,1,1.5,2,2.5,3,3.75,4,5])
w=2*np.pi/T
sp = [s for s in os.listdir(ss) if (s.startswith('profile'))]
gms = [s for s in os.listdir(ss) if (s.endswith('.AT2'))]
# gms = [s for s in os.listdir(ss) if (s.startswith('RSN'))]
gms = sorted(gms)
ii6 = 0
iii = 0

SF = np.arange(0.25,10.5,0.25)
# np.arange(0.25,10.5,0.25)

for ii in np.arange(0,len(SF),1):
    SF_req = SF[ii]
    ii6 = 0

    #result = []
    #SA = []
    #SA2 = []
    #gmotion = []
    #RT = []
    for s in sp:
        print(s)
        sprofile = EquivLin.input_soil_profile(s)  # Better way
        for g in gms:
            iii = iii + 1
            gmotion = EquivLin.GroundMotion(filename=g, Butrwrth=True, scalefactor = SF_req)
            result = EquivLin.Dynamics(sprofile,gmotion, Error=error, modtype='FreqInd',verbose=False, run_all=False, strainratio=0.65, iters = 30)# 0.65
            print(result.Dmin)
            if result.ReportError is 'n':
                taumax,gammax = result._calcTauGam(MaxOnly=True,domain='time',returns=True)
                SA, w = result.RespSpec(0,w = w,MotType='outcrop')
                SA2, w = result.RespSpec(sum(sprofile.t)+0.1,w = w,MotType='outcrop')

                if ii6 is 0:
                    g_max= max(gammax)
                    RT = SA/SA2
                    # RT
                    SAr=SA2;
                    form = '%16.12f'
                else:
                    R = SA/SA2
                    # R
                    RT = np.column_stack((RT, R))
                    SAr = np.column_stack((SAr, SA2))
                    g_max = np.column_stack((g_max, max(gammax)))
                    form += '%16.12f'
                ii6+=1
            else:
                neg = iii

    f = open('AmpFac_'+str(SF_req)+'.txt', 'w')
    np.savetxt(f, RT, fmt=form, delimiter=' ', newline='\n', header='', footer='', comments='# ')
    f.close()

    f = open('SaRock_'+str(SF_req)+'.txt', 'w')
    np.savetxt(f, SAr, fmt=form, delimiter=' ', newline='\n', header='', footer='', comments='# ')
    f.close()

    print('Progress:' + str((ii+1)/len(SF)))

# f = open('MOTION_Vsgam.txt', 'w')
# np.savetxt(f, g_max, fmt=form, delimiter=' ', newline='\n', header='', footer='', comments='# ')
# f.close()
