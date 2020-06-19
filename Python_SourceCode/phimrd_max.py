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
# T=np.array([0.01,0.0108,0.0117,0.0126,0.0136,0.0147,0.0158,0.0171,0.0185,0.02,0.0215,0.0233,0.0251,0.0271,0.0293,0.0316,0.0341,0.0369,0.0398,0.0430,0.0464,0.0501,0.0541,0.0584,0.0631,0.0681,0.0736,0.0794,0.0858,0.0926,0.1,0.1080,0.1166,0.1259,0.1359,0.1468,0.1585,0.1711,0.1848,0.1955,0.2154,0.2326,0.2512,0.2712,0.2929,0.3162,0.3415,0.3687,0.3981,0.4299,0.4642,0.5012,0.5412,0.5843,0.6310,0.6813,0.7356,0.7943,0.8577,0.9261,1,1.0798,1.1659,1.2589,1.3594,1.4678,1.5849,1.7113,1.8478,1.9953,2.1544,2.3263,2.5119,2.7123,2.9286,3.1623,3.4145,3.6869,3.9811,4.2987,4.6416,5.0119,5.4117,5.8434,6.3096,6.8129,7.3564,7.9433,8.5770,9.2612,10])
w=2*np.pi/T
sp = [s for s in os.listdir(ss) if (s.startswith('profile'))]
gms = [s for s in os.listdir(ss) if (s.endswith('.AT2'))]
# gms = [s for s in os.listdir(ss) if (s.startswith('RSN'))]
gms = sorted(gms)
ii6 = 0
iii = 0
for s in sp:
    print(s)
    sprofile = EquivLin.input_soil_profile(s)  # Better way
    for g in gms:
        iii = iii + 1
        gmotion = EquivLin.GroundMotion(filename=g, Butrwrth=True, scalefactor = 1)
        result = EquivLin.Dynamics(sprofile,gmotion, Error=error, modtype='FreqInd',verbose=False, run_all=True, strainratio=0.65)# 0.65
        print(result.Dmin)
        if result.ReportError is 'n':
            taumax,gammax = result._calcTauGam(MaxOnly=True,domain='time',returns=True)
            SA, w = result.RespSpec(0,w = w,MotType='outcrop')
            SA2, w = result.RespSpec(sum(sprofile.t)+0.1,w = w,MotType='outcrop')
    
            if ii6 is 0:
               g_max= max(gammax)
               RT = SA/SA2
               RT
               SAr=SA2;
               form = '%16.12f'
            else:
                R = SA/SA2
                R
                RT = np.column_stack((RT, R))
                SAr = np.column_stack((SAr, SA2))
                g_max = np.column_stack((g_max, max(gammax)))
                form += '%16.12f'
            ii6+=1
        else:
            neg = iii    

f = open('Amp_Fac.txt', 'w')
np.savetxt(f, RT, fmt=form, delimiter=' ', newline='\n', header='', footer='', comments='# ')
f.close()

f = open('SA_rock.txt', 'w')
np.savetxt(f, SAr, fmt=form, delimiter=' ', newline='\n', header='', footer='', comments='# ')
f.close()

# f = open('MOTION_Vsgam.txt', 'w')
# np.savetxt(f, g_max, fmt=form, delimiter=' ', newline='\n', header='', footer='', comments='# ')
# f.close()
