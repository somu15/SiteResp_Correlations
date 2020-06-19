
#!/usr/bin/python
# IZ_1993.py
from __future__ import division
'''
    Copyright (C) 2014 Samuel J. Lasley, Russell A. Green, and
						Adrian Rodriguez-Marek

    Permission is hereby granted, free of charge, to any person 
    obtaining a copy of this software and associated documentation 
    files (the "Software"), to deal in the Software without restriction, 
    including without limitation the rights to use, copy, modify, merge, 
    publish, distribute, sublicense, and/or sell copies of the Software, 
    and to permit persons to whom the Software is furnished to do so, 
    subject to the following conditions:

    The above copyright notice and this permission notice shall be 
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
    SOFTWARE.
'''
#Edited at 9/14/2016
#import pdb
import numpy as np
import matplotlib.pyplot as plt
#straindef=np.logspace(-6,-0.004364805,num=10)

straindef = np.array([0.000001, 0.000003, 0.00001, 0.00003, 0.0001, 0.0003,
                      0.001, 0.003, 0.01, 0.99, 5])


def IZ_1993(gam=straindef, PI=0., sigm=100):
    ''' Returns the Ishibashi and Zhang (1993) values for G/Gmax
    and damping (decimal) for a given values of strain (as a decimal),
    Plasticity index, and initial mean effective confining stress (kPa).
    '''

    PI = float(PI)
    sigm = float(sigm)

    if PI == 0.:
        n = 0.0
    elif PI <= 15:
        n = 3.37e-6 * PI ** 1.404
    elif PI <= 70:
        n = 7e-7 * PI ** 1.976
    else:
        n = 2.7e-5 * PI ** 1.115


    K = 0.5 * (1 + np.tanh(np.log(((0.000102 + n) / gam) ** 0.492)))
    m = (0.272 * (1 - np.tanh(np.log((0.000556 / gam) ** 0.4))) *
                np.exp(-0.0145 * PI ** 1.3))
    Gratio = K * sigm ** m
    np.putmask(Gratio, Gratio>1, 1)
    #use np.copyto when updates to numpy 1.7 or use Gratio.clip(max=1)
    D = 0.1665 * (1 + np.exp(-0.0145 * PI ** 1.3)) * (0.586 * Gratio ** 2. - 
                1.547 * Gratio + 1.)
    Dmin=0.1665 * (1 + np.exp(-0.0145 * PI ** 1.3)) * (0.586 * 1 ** 2. - 
                1.547 * 1 + 1.)
    return (gam, Gratio, D, Dmin)


def DS_2001(gam=straindef, PI=0., sigm=100., OCR=1., soil=0, N=10., frq=1.,MRD_rand=0,epsilon=-10,Su=-1.0,Gmax=375000):
    ''' 
    Darendeli and Stokoe (2001) shear modulus reduction and damping curves.
    gam= shear strain (decimal)
    PI= Plasticity Index
    sigm = the initial mean effective confining stress (kPa)
    OCR = the overconsolidation ratio
    soil = the soil type indicator:
        0 for general
        1 for clean sand
        2 for sands with high fines contents
        3 for silts
        4 for clays
    #Mahdi's change
    MRD_rand = Randomizing MRD curves
        0 Use median MRD curve
        1 Generate a randmized MRD curve
    epsilon = The random number used in generating MRD realizations
        -10 The code produces a random number
        Input a two dimentional vector to specify it manually
    Su = shear strength of soil( To not correct NG for shear strenght input -1)
    Gmax = Maximum soil modulus
    #
    N = the number of cycles
    frq = the frequency of loading
    Dmintype = 1 for 'Default', 0 for 'Green'  (Suggest using Green's for now!!!!!)
                !!!Not Implemented !!!
    The output of the function is:
    gam = the same input shear strain as a decimal
    Gratio = G/Gmax = the shear modulus normalized by the small strain
        shear modulus.
    D = damping (decimal)
    '''
    PI = float(PI)
    sigm = float(sigm)
    OCR = float(OCR)
    N = float(N)
    frq = float(frq)
    soil = int(soil)
    MRD_rand = int(MRD_rand)
    Su=float(Su)
    Gmax=float(Gmax)
    phi1 = np.array([ 3.52E-02, 4.74E-02, 3.34E-02, 4.16E-02, 2.58E-02])
    phi2 = np.array([ 1.01E-03,-2.34E-03,-5.79E-05, 6.89E-04, 1.95E-03])
    phi3 = np.array([ 3.246E-01, 2.50E-01, 2.49E-01, 3.21E-01, 9.92E-02])
    phi4 = np.array([ 3.483E-01, 2.34E-01, 4.82E-01, 2.80E-01, 2.26E-01])
    phi5 = np.array([ 9.19E-01, 8.95E-01, 8.45E-01, 1.00E+00, 9.75E-01])
    phi6 = np.array([ 8.005E-01, 6.88E-01, 8.89E-01, 7.12E-01, 9.58E-01])
    phi7 = np.array([ 1.29E-02, 1.22E-02, 2.02E-02, 3.03E-03, 5.65E-03])
    phi8 = np.array([-1.069E-01, -1.00E-01,-1.00E-01,-1.00E-01,-1.00E-01])
    phi9 = np.array([-2.889E-01, -1.27E-01,-3.72E-01,-1.89E-01,-1.96E-01])
    phi10 = np.array([ 2.919E-01, 2.88E-01, 2.33E-01, 2.34E-01, 3.68E-01])
    phi11 = np.array([ 6.329E-01, 7.67E-01, 7.76E-01, 5.92E-01, 4.66E-01])
    phi12 = np.array([-5.66E-03, -2.83E-02,-2.94E-02,-7.67E-04, 2.23E-02])
    phi13 = -4.23
    phi14 = 3.62 
    phi15 = -5.00
    phi16 = -2.50E-01
    
    # Calculate the reference shear strain:
    gamr = ((phi1[soil] + phi2[soil] * PI * OCR ** phi3[soil]) * (sigm / 101.325) ** phi4[soil])
    gam = gam * 100  # convert to percent
    
    # Calculate assorted variables:
    a = phi5[soil]
    b = phi11[soil] + phi12[soil] * np.log(N)
    C1 = -1.1143 * a ** 2 + 1.8618 * a + 0.2523
    C2 = 0.0805 * a ** 2 - 0.0710 * a - 0.0095
    C3 = -0.0005 * a ** 2 + 0.0002 * a + 0.0003
    # Calculate G/Gmax
    Gratio = 1. / (1. + (gam / gamr) ** a)

    
    Dmasinga1 = 100 / np.pi * (4 * (gam - gamr * np.log((gam + gamr) 
                / gamr)) / (gam ** 2 / (gam + gamr)) - 2 )
                
    Dmasing = C1 * Dmasinga1 + C2 * Dmasinga1 ** 2 + C3 * Dmasinga1 ** 3
    
    Dmin = ((phi6[soil] + phi7[soil] * PI * OCR ** phi8[soil]) * (sigm / 101.3) ** phi9[soil] * (1 + phi10[soil] * np.log(frq)))
    
    D = b * (Gratio ** 0.1) * Dmasing + Dmin


    
    if MRD_rand == 1:
        ro=-0.5 # <--- This value is based on a sentence in darandali thesis
        sNG = np.exp(phi13)+np.sqrt(0.25/np.exp(phi14)-((Gratio-0.5)**2)/np.exp(phi14))
        sD = np.exp(phi15)+np.exp(phi16)*np.sqrt(D)
        sDmin=np.exp(phi15)+np.exp(phi16)*np.sqrt(Dmin)
        # Add varibility due to PI
        dDmas_dgr = 8.74/(gamr/100);
        PI_v=15 #<----
        Dmin_v= (phi7[soil]*OCR**phi7[soil]*(sigm/101)**(phi9[soil])*(1+phi10[soil]*np.log(frq)))*PI_v;
        gr_v=(phi2[soil]*OCR**phi3[soil])*(sigm/101)**phi4[soil]/100*PI_v;
        NG_gr=a/(4*(gamr/100))*gr_v;# Derivative of NG with respect to gr
        D2_v = abs(b*0.1*1.866*(-a/(4*(gamr/100)))*13.57+b*(1/2)**0.1*dDmas_dgr)*gr_v
        sD_Dmin =np.exp(phi15)+np.exp(phi16)*np.sqrt(Dmin);
        s2D2_gr = np.exp(phi16)**2*(b*(0.5)**0.1*13.57)+2*np.exp(phi15)*np.exp(phi16)*(np.sqrt(b*(0.5)**0.1*13.57+Dmin)-np.sqrt(Dmin))
        s2D = sD_Dmin**2+Dmin_v**2+((s2D2_gr+D2_v)/s2D2_gr)*(np.exp(phi16)**2*(D-Dmin)+2*np.exp(phi15)*np.exp(phi16)*(np.sqrt(D)-np.sqrt(Dmin)))
        sD=np.sqrt(s2D)
        sNG = np.exp(phi13)+np.sqrt(4*NG_gr**2+1/np.exp(phi14))*np.sqrt(0.25-((Gratio-0.5)**2))
        # Make the variability of G/Gmax before threshold strain zero 
        gtf_r=np.array([0.000331314,0.000479492,0.000700539,0.00095242,0.001268053,0.001704959,0.00220124,0.002564315,0.00389172,0.00389172])
        PI_r=np.array([0,4.637925427,11.70840849,18.02339908,25.00172487,33.18696591,39.62171166,44.97718245,60,200])
        gtf=np.exp(np.interp(PI,PI_r,np.log(gtf_r)))
        ga=gtf*30
        z=np.ones(len(gam))
        for i in range(len(gam)):
            if gam[i]>=ga:
                z[i]=1
            elif gam[i]<gtf:
                z[i]=0
            else:
                z[i] = ((np.log(gam[i])-np.log(gtf))/(np.log(ga)-np.log(gtf)))
                
        sNG=z*sNG
        # Bring the variability and sigma from lognormal to normal distribution
        D_log = np.log(D/np.sqrt(1+sD**2/D**2))
        Dmin_log = np.log(Dmin/np.sqrt(1+sDmin**2/Dmin**2))
        Gratio_log = np.log(Gratio/np.sqrt(1+sNG**2/Gratio**2))
        sD_log = np.sqrt(np.log(1+sD**2/D**2))
        sDmin_log = np.sqrt(np.log(1+sDmin**2/Dmin**2))
        sNG_log = np.sqrt(np.log(1+sNG**2/Gratio**2))
        
        if epsilon==-10: #If no epsilon is assigned, make it. If it is assigned, just read them
            epsilon1=-5
            epsilon2=-5
            while abs(epsilon1)>=2 or abs(epsilon1)>=2:
                epsilon1 = np.random.normal(loc=0.0, scale=1.0)
                epsilon2 = np.random.normal(loc=0.0, scale=1.0)
                epsilon2=((1-ro**2)**0.5)*epsilon2+ro*epsilon1
        else:
            epsilon1=epsilon[0]
            epsilon2=epsilon[1]
        Gratio_n_log=Gratio_log+epsilon1*sNG_log
        D_n=np.exp(D_log+epsilon2*sD_log)
        Dmin_n=np.exp(Dmin_log+epsilon2*sDmin_log)
        
                
        D=D_n
        Dmin=Dmin_n
        D=np.minimum(D,15)
        Gratio = np.exp(Gratio_n_log)


        # If you have Su value correct the curve for it
        if Su!=-1.0:
            sSu = 0.5 #<----
            
            # Yee et al. (2013)
            g1_vector = [0.3, 0.1, 0.07, 0.05, 0.03, 0.01, 0.005,0.001]
            
            for g1 in g1_vector:
                    
                
                condition = gam > g1

                A = np.extract(condition , gam)

                g_1 = A[0]

                A = np.extract(condition , Gratio)

                Gratio_1 = A[0]

                condition = gam<g1

                A = np.extract(condition , gam)

                g_2 = A[-1]

                A = np.extract(condition ,Gratio )

                Gratio_2 = A[-1]

                NG1=(Gratio_2*g_2 - Gratio_1*g_1) / (g_2 - g_1)

                tauf_n=np.exp(np.log(Su) + epsilon1 *sSu )

                gref=(tauf_n - (g_1/100)*Gratio_1*Gmax)/(NG1*Gmax)*100

                if gref>0 and NG1>0:
                    break

            count=0

            for g in gam:

                if g>g1:
                    
                    Gratio[count] = (g_1*Gratio_1+np.sign(gref)*NG1*(g-g_1)/(1+(g-g_1)/np.absolute(gref)))/g

                count+=1

                    




    elif Su != -1.0:
        # Yee et al. (2013)
        g1_vector = [0.3, 0.1, 0.07, 0.05, 0.03, 0.01, 0.005]

        for g1 in g1_vector:

            condition = gam > g1

            A = np.extract(condition , gam)

            g_1 = A[0]

            A = np.extract(condition , Gratio)

            Gratio_1 = A[0]

            condition = gam<g1

            A = np.extract(condition , gam)

            g_2 = A[-1]

            A = np.extract(condition ,Gratio )

            Gratio_2 = A[-1]

            NG1=(Gratio_2*g_2 - Gratio_1*g_1) / (g_2 - g_1)

            tauf_n=Su

            gref=(tauf_n - (g_1/100)*Gratio_1*Gmax)/(NG1*Gmax)*100
            
            if gref>0 and NG1>0:
                break
        
        count=0


        for g in gam:

            if g>g1:
                
                Gratio[count] = (g_1*Gratio_1+np.sign(gref)*NG1*(g-g_1)/(1+(g-g_1)/np.absolute(gref)))/g


            count+=1

    return gam / 100, Gratio, D/100 , Dmin/100
    

def Menq(gam = straindef, Cu = 7 , D50=1 , sigm=100.,N = 10., frq = 1.,MRD_rand=0,epsilon=-10,Su=-1,Gmax=375000):
    '''    
    Usage: [NG,D ,sNG, sD] = Menq(g,so
    g = strain vector (enter a negative number to use default)
    unitless (e.g., not in percent)
    Cu   Coeff. of Uniformity
    D50 mean grain size diam (in mm)
    so  mean effective stress(in kPa)
    #Mahdi's change
    MRD_rand = Randomizing MRD curves
        0 Use median MRD curve
        1 Generate a randmized MRD curve
    epsilon = The random number used in generating MRD realizations
        -10 The code produces a random number
        Input a two dimentional vector to specify it manually
    Su = shear strength of soil( To not correct NG for shear strenght input -1)
    Gmax = Maximum soil modulus
    #
    The last two argument are optional (default is used otherwise)
    Ncyc: number of cycles of loading (default is 10)
    f: frequency (Hz) (default is 1 Hz)
    Output
    NG: normalized modulus reduction
    D: Damping curve
    sNG: Standard deviation of normalized curve
    sD: Standard deviation of damping curve
    '''


    gam=gam*100;
    pa=101.325;
    gamr = 0.12 * Cu ** (-0.6) * (sigm/pa) ** (0.5*Cu ** (-0.15));
    a = 0.86 + 0.1*np.log10(sigm/pa);
    Gratio=1/(1+(gam/gamr)**a);

    b=0.6329-0.0057*np.log(N);
    


    c1=-1.1143*a**2 + 1.8618*a + 0.2523;
    c2=0.0805* a**2 - 0.0710*a - 0.0095;
    c3=-0.0005*a**2 + 0.0002*a + 0.0003;
    Dmas_a1=(100/3.14*(4*((gam-gamr*np.log((gam+gamr)/gamr))/(gam**2/(gam+gamr)))-2));
    Dmas = c1*Dmas_a1 + c2*Dmas_a1**2 + c3*Dmas_a1**3;
    Dmin = 0.55*Cu**0.1*D50**(-0.3)*(sigm/pa)**(-0.05);

    #Dec 18,2015: change the coefficient of so/pa to -0.08 from -.05. this is 
    #a discrepancy in Menq's dissertation
    #Dec 21, 2015: email from Menq
    #Dear Adrian,
    #Happy Holidays!
    # The value of -0.08 was the best fit for my data.  However, I was concerned that I might be over fitting the data. So, I rounded up the number to -0.05, and it worked just fine. The value of -0,05 was what I recommended.  I apologized for the writing.  It was painful for me to write in English. 
    # Please note that most of the test specimens are air dried.   With added moisture, Dmin is expected to be higher. 
    # Best regards,
    #
    #
    # Menq, 


    F = b* Gratio**0.1
    D = F*Dmas + Dmin
    phi13 = -4.23
    phi14 = 3.62 
    phi15 = -5.00
    phi16 = -2.50E-01
    if MRD_rand == 1:
        ro=-0.5
        sNG = np.exp(phi13)+np.sqrt(0.25/np.exp(phi14)-((Gratio-0.5)**2)/np.exp(phi14))
        sD = np.exp(phi15)+np.exp(phi16)*np.sqrt(D)
        sDmin=np.exp(phi15)+np.exp(phi16)*np.sqrt(Dmin)
        PI=0
        gtf_r=np.array([0.000331314,0.000479492,0.000700539,0.00095242,0.001268053,0.001704959,0.00220124,0.002564315,0.00389172,0.00389172])
        PI_r=np.array([0,4.637925427,11.70840849,18.02339908,25.00172487,33.18696591,39.62171166,44.97718245,60,200])
        gtf=np.exp(np.interp(PI,PI_r,np.log(gtf_r)))
        ga=gtf*30
        z=np.ones(len(gam))
        for i in range(len(gam)):
            if gam[i]>=ga:
                z[i]=1
            elif gam[i]<gtf:
                z[i]=0
            else:
                z[i] = ((np.log(gam[i])-np.log(gtf))/(np.log(ga)-np.log(gtf)))
        sNG=z*sNG
                
        D_log = np.log(D/np.sqrt(1+sD**2/D**2))
        Dmin_log = np.log(Dmin/np.sqrt(1+sDmin**2/Dmin**2))
        Gratio_log = np.log(Gratio/np.sqrt(1+sNG**2/Gratio**2))
        sD_log = np.sqrt(np.log(1+sD**2/D**2))
        sDmin_log = np.sqrt(np.log(1+sDmin**2/Dmin**2))
        sNG_log = np.sqrt(np.log(1+sNG**2/Gratio**2))
        if epsilon==-10:
            epsilon1=-5
            epsilon2=-5
            while abs(epsilon1)>=2 or abs(epsilon1)>=2:
                epsilon1 = np.random.normal(loc=0.0, scale=1.0)
                epsilon2 = np.random.normal(loc=0.0, scale=1.0)
                epsilon2=((1-ro**2)**0.5)*epsilon2+ro*epsilon1
        else:
            epsilon1=epsilon[0]
            epsilon2=epsilon[1]
        Gratio_n=np.exp(Gratio_log+epsilon1*sNG_log)
        D_n=np.exp(D_log+epsilon2*sD_log)
        Dmin_n=np.exp(Dmin_log+epsilon2*sDmin_log)    
                
        D=D_n
        Gratio=Gratio_n
        Dmin=Dmin_n
        D=np.minimum(D,15)

    
    return gam / 100, Gratio, D/100 , Dmin/100

def PietPeat(gam= straindef, sigm = 100 ,N = 10., frq = 1.,MRD_rand=0,epsilon=-10,Su=-1,Gmax=375000):
    '''
       Piet fit to Peat data
       Usage: [NG,D ,sNG, sD] = Darendeli(g,PI,svo,Ncyc,f)
       Input
           g = strain vector (enter a negative number to use default)
               unitless (e.g., not in percent)
           PI?????
           so  mean effective stress(in kPa)
           OCR
           Ko: lateral earth pressure coeff. at rest
         The last two argument are optional (default is used otherwise)
           Ncyc: number of cycles of loading (default is 10)
           f: frequency (Hz) (default is 1 Hz)
       Output
           NG: normalized modulus reduction
           D: Damping curve
           sNG: Standard deviation of normalized curve
           sD: Standard deviation of damping curve
    '''

    #enter parameters of model, from Table 8.12 of Darendeli
    p=np.array([3.52E-02,1.01E-03,3.25E-01,3.48E-01,9.19E-01,8.01E-01 ,1.29E-02,-1.07E-01,-2.89E-01,2.92E-01,6.33E-01,-5.66E-03, -4.23,3.62,-5.00,-2.50E-01,5.62,2.78]);
    gam=gam

    pa = 101.3;
    
    a=0.776;

    gamr=(0.995*(sigm/pa)**0.694)/100;

    c1=-1.1143*a**2 + 1.8618*a + 0.2523;
    c2=0.0805* a**2 - 0.0710*a - 0.0095;
    c3=-0.0005*a**2 + 0.0002*a + 0.0003;

    Gratio=1/(1+(gam/gamr)**a);

    Dmas_a1=(100/3.14*(4*((gam-gamr*np.log((gam+gamr)/gamr))/(gam**2./(gam+gamr)))-2));
    Dmas = c1*Dmas_a1 + c2*Dmas_a1**2 + c3*Dmas_a1**3;
    Dmin =  2.512*(sigm/pa)**(-0.2889);
    F = 0.712* Gratio**0.1;
    # F=0.637;
    D = F*Dmas + Dmin
    phi13 = -4.23
    phi14 = 3.62 
    phi15 = -5.00
    phi16 = -2.50E-01
    if MRD_rand == 1:
        PI=50
        ro=-0.5
        sNG = np.exp(phi13)+np.sqrt(0.25/np.exp(phi14)-((Gratio-0.5)**2)/np.exp(phi14))
        sD = np.exp(phi15)+np.exp(phi16)*np.sqrt(D)
        sDmin=np.exp(phi15)+np.exp(phi16)*np.sqrt(Dmin)
        # Make the variability of G/Gmax before threshold strain zero 
        gtf_r=np.array([0.000331314,0.000479492,0.000700539,0.00095242,0.001268053,0.001704959,0.00220124,0.002564315,0.00389172,0.00389172])
        PI_r=np.array([0,4.637925427,11.70840849,18.02339908,25.00172487,33.18696591,39.62171166,44.97718245,60,200])
        gtf=np.exp(np.interp(PI,PI_r,np.log(gtf_r)))
        ga=gtf*30
        z=np.ones(len(gam))
        for i in range(len(gam)):
            if gam[i]*100>=ga:
                z[i]=1
            elif gam[i]*100<gtf:
                z[i]=0
            else:
                z[i] = ((np.log(gam[i]*100)-np.log(gtf))/(np.log(ga)-np.log(gtf)))
        sNG=z*sNG
                
        D_log = np.log(D/np.sqrt(1+sD**2/D**2))
        Dmin_log = np.log(Dmin/np.sqrt(1+sDmin**2/Dmin**2))
        Gratio_log = np.log(Gratio/np.sqrt(1+sNG**2/Gratio**2))
        sD_log = np.sqrt(np.log(1+sD**2/D**2))
        sDmin_log = np.sqrt(np.log(1+sDmin**2/Dmin**2))
        sNG_log = np.sqrt(np.log(1+sNG**2/Gratio**2))
        if epsilon==-10:
            epsilon1=-5
            epsilon2=-5
            while abs(epsilon1)>=2 or abs(epsilon1)>=2:
                epsilon1 = np.random.normal(loc=0.0, scale=1.0)
                epsilon2 = np.random.normal(loc=0.0, scale=1.0)
                epsilon2=((1-ro**2)**0.5)*epsilon2+ro*epsilon1
        else:
            epsilon1=epsilon[0]
            epsilon2=epsilon[1]
        Gratio_n=np.exp(Gratio_log+epsilon1*sNG_log)
        D_n=np.exp(D_log+epsilon2*sD_log)
        Dmin_n=np.exp(Dmin_log+epsilon2*sDmin_log)    
                
        D=D_n
        Gratio=Gratio_n
        Dmin=Dmin_n
        D=np.minimum(D,15)

        # If you have Su value correct the curve for it
        if Su!=-1.0:
            
            sSu = 0.5 #<---
            # Yee et al. (2013)
            g1_vector = np.array([0.3, 0.1, 0.07, 0.05, 0.03, 0.01, 0.005]) / 100
            for g1 in g1_vector:
                
                condition = gam > g1

                A = np.extract(condition , gam)

                g_1 = A[0]

                A = np.extract(condition , Gratio)

                Gratio_1 = A[0]

                condition = gam<g1

                A = np.extract(condition , gam)

                g_2 = A[-1]

                A = np.extract(condition ,Gratio )

                Gratio_2 = A[-1]

                NG1=(Gratio_2*g_2 - Gratio_1*g_1) / (g_2 - g_1)

                tauf_n=np.exp(np.log(Su) + epsilon1 *sSu )

                gref=(tauf_n - (g_1/100)*Gratio_1*Gmax)/(NG1*Gmax)*100

                if gref>0 and NG1>0:
                    
                    break
            


            count=0

            for g in gam:

                if g>g1:
                    
                    Gratio[count] = (g_1*Gratio_1+np.sign(gref)*NG1*(g-g_1)/(1+(g-g_1)/np.absolute(gref)))/g

                count+=1

                    




    elif Su != -1.0:

        # Yee et al. (2013)
        g1_vector =np.array([0.3, 0.1, 0.07, 0.05, 0.03, 0.01, 0.005]) / 100

        for g1 in g1_vector:

            condition = gam > g1

            A = np.extract(condition , gam)

            g_1 = A[0]

            A = np.extract(condition , Gratio)

            Gratio_1 = A[0]

            condition = gam<g1

            A = np.extract(condition , gam)

            g_2 = A[-1]

            A = np.extract(condition ,Gratio )

            Gratio_2 = A[-1]

            NG1=(Gratio_2*g_2 - Gratio_1*g_1) / (g_2 - g_1)

            tauf_n=Su

            gref=(tauf_n - (g_1/100)*Gratio_1*Gmax)/(NG1*Gmax)*100
            
            if gref>0 and NG1>0:
                break

        count=0


        for g in gam:

            if g>g1:
                Gratio[count] = (g_1*Gratio_1+np.sign(gref)*NG1*(g-g_1)/(1+(g-g_1)/np.absolute(gref)))/g


            count+=1            
    return gam , Gratio, D/100 , Dmin/100




def Vucetic_Dobry (gam = straindef, PI=0):
    '''
    PI = Plasticity Index
    '''
    gam=gam*100
    if PI==0:
        gamref = np.array ([1.0e-4, 3.16e-4 ,1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1])
        Gref = np.array ([1 , 1 , 0.96 , 0.88 , 0.7 , 0.47 , 0.26 , 0.11 , 0.03])
        Dref = np.array ([1 , 1 , 1 , 3 , 5.4 , 9.8 , 15 , 20.3 , 24])
    elif PI==15:
        gamref = np.array ([1.0e-4, 3.16e-4 ,1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1])
        Gref = np.array ([1 , 1 , 0.99 , 0.94 , 0.81 , 0.64 , 0.41 , 0.22 , 0.1])
        Dref = np.array ([1 , 1 , 1 , 2.6 , 4.5 , 7.5 , 11.6 , 16 , 20])
    elif PI==30:
        gamref = np.array ([1.0e-4, 3.16e-4 ,1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1])
        Gref = np.array ([1 , 1 , 1 , 0.98 , 0.9 , 0.75 , 0.53 , 0.35 , 0.17])
        Dref = np.array ([1 , 1 , 1 , 2.1 , 3.8 , 5.9 , 8.8 , 12.5 , 16.9])
    elif PI==50:
        gamref = np.array ([1.0e-4, 3.16e-4 ,1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1])
        Gref = np.array ([1 , 1 , 1 , 1 , 0.95 , 0.84 , 0.67 , 0.47 , 0.25])
        Dref = np.array ([1 , 1 , 1 , 1.8 , 2.9 , 4.3 , 6.2 , 9.5 , 13.5])
    elif PI==100:
        gamref = np.array ([1.0e-4, 3.16e-4 ,1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1])
        Gref = np.array ([1 , 1 , 1 , 1 , 0.98 , 0.92 , 0.81 , 0.63 , 0.37])
        Dref = np.array ([1 , 1 , 1 , 1.5 , 2 , 2.9 , 4.1 , 6.2 , 9.8])
    elif PI==200:
        gamref = np.array ([1.0e-4, 3.16e-4 ,1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1])
        Gref = np.array ([1 , 1 , 1 , 1 , 1 , 0.96 , 0.89 , 0.75 , 0.53 ])
        Dref = np.array ([1 , 1 , 1.1 , 1.3 , 1.6 , 2.1 , 3 , 4.8 , 8.1 ])
    else:
        print('In Vucetic_Dobry model PI should be one of these values: 0,15, 30,50,100,200')
    Gratio = np.interp( np.log10(gam) , np.log10(gamref), Gref)
    D = np.interp( np.log10(gam) , np.log10(gamref), Dref)
    return gam / 100, Gratio, D/100




def Undamped_linear (gam = straindef):
    '''
    PI = Plasticity Index
    '''

    Gratio = np.ones(len(gam))
    D = np.zeros(len(gam))
    print(Gratio)
    return gam , Gratio, D/100 

