#!/usr/bin/python
# -*- coding: utf-8 -*-
# EquivLin.py

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

# Edited at 9/14/2016
import csv
import pdb
import numpy as np
#import numexpr as ne
from scipy.interpolate import interp1d
#import guiqwt.pyplot as plt
import matplotlib.pyplot as plt
import degrcurv as dc

#import wxfiledialog as fd
#import qtfiledialog as fd
#import cPickle
try:
    xrange
except NameError:
    xrange = range

class SoilProfile(object):
    '''
    This object acts as a container for the Soil Profile and subsequent
        calculations.

    Example:

    SoilProfileObject = SoilProfile(profile_name, Unit_Wt, Vs, t, **kwargs)

    The functions 'inputread' and 'input_soil_profile' automatically call this
    object.  Alternately, the required values can be given manually.

        Inputs:
            profile_name = a name for the profile
            Unit_wt = an array of unit weights of the layers of the profile.
            Vs = an array of shear wave velocities of the layers of the profile.
            t = tan array of the thicknesses of each layer of the profile.

        Keyword Arguments:           Definition/Default:
            filename         filename of the input profile

            ko               The ratio of vertical over horizontal shear
                                stresses (array). This is used to calculate
                                the degradation curves.

            phi              The friction angle (array). Also used to calculate
                                degradation curves. Used if the ko array
                                is zero.

            N160             Corrected SPT blowcounts (array).  Only used if ko
                                and phi are zero. If an N160 value is not
                                given, ko is assumed to be 0.5.

            gwt              The depth of the ground water table.  If not given,
                                it is assumed to be at the surface.

            g                Acceleration due to gravity.
                                Assumed to be 9.81 m/s^2.

            layerno          An array of integers numbering the layers.
                                Assumed to be counted from the surface
                                starting at 1.

           layername         An array of names describing the layers.
           scalefactor       A factor by which ground motions can be amplified
                               or deamplified. Assumed to be 1.

           strainratio      Ratio of equivalent uniform strain divided by the
                               maximum strain. Default is 0.65. This value is
                               the default when calling the
                               Dynamics.iterations method.

           loc_mot          Layer at which the motion is input into the
                               profile, given as an integer.  When using this
                               value, the Dynamics._calcWaves method counts the
                               layers from the surface, starting at 1,
                               irrespective of the layerno array values.
                               Default is the bottom layer (bedrock).

           outcrop_mot      0 signifies that the input motion is an outcrop
                               motion, 1 otherwise. This can be overridden by
                               the Dynamics object call if the 'outcrop_mot'
                               keyword is used.
           OCR              An array giving the Over Consolidation Ratio for
                               each layer.  This is used only for the
                               Darendeli and Stokoe (2001) Degradation Curves.

           DS_soil          Soil number for Darendeli and Stokoe Degrad Curves
                               The soil type indicator:
                                    0 for a general set of curves
                                     (Table 8.12 from Darendeli's Dissertation)
                                    1 for clean sand
                                    2 for sands with high fines contents
                                    3 for silts
                                    4 for clays

           DS_N             Number of cycles for the Darendeli and
                               Stokoe (2001) Degradation Curves.
                               Default is 10 for all layers.

           DS_freq          Frequency for the Darendeli and
                               Stokoe (2001) Degradation Curves.
                               Default is 1 Hz for all layers.
 ======================   NOT IMPLEMENTED =====================================
           ACCTHo6          an array indicating if the acceleration time history should be
                               output for the center of each layer.  1 will output, 0 will suppress output.
           AccTHoutcrop6    Specifies if the motion from the above is an outcrop motion (O) or within (1)[default]
           AccTHotype6      Specifies if the above timehistories should be the entire time history (1),
                              or the max value only (0)[default]
          TauTH7a           An array specifying if the stress time histories of the corresponding layers should be output.
                              1 will output, 0 will suppress output.
          GamTH7b           An array specifying if the stress time histories of the corresponding layers should be output.
                               1 will output, 0 will suppress output.
          RespSpec9         An array specifying if the response spectrum should be calculated and output:
                              0 for no [default], 1 for outcrop, 2 for within


          There should be additional options that I need to code.<<<+========================================================


    The following 'data attributes' or variables are created by __init__:
        self.t
        self.unit_wt
        self.gwt
        self.profile_name
        self.ko
        self.phi
        self.PI
        self.N160
        self.g
        self.layerno
        self.layername
        self.scalefactor
        self.strainratio
        self.loc_mot
        self.outcrop_mot
        self.filename
        self.OCR

        Some of these may equal 'None'.
    '''
    def __init__(self, profile_name, unit_wt, Vs, t, **kwargs):
        self.profile_name = profile_name
        self.unit_wt = unit_wt
        self.Vs = Vs
        self.t = t
        self.units = kwargs.get('units',0)
        self.ko = kwargs.get('ko', np.zeros(len(t),))
        self.phi = kwargs.get('phi', np.zeros(len(t),))
        self.PI = kwargs.get('PI', np.zeros(len(t),))
        self.N160 = kwargs.get('N160', np.zeros(len(t),))
        self.gwt = kwargs.get('gwt', 0)
        self.g = kwargs.get('g', 9.81)
        self.layerno = kwargs.get('layerno', np.arange(len(t),)+1)
        self.layername = kwargs.get('layername')
        self.scalefactor = kwargs.get('scalefactor',1.)
        self.strainratio = kwargs.get('strainratio',0.65)
        self.loc_mot = kwargs.get('loc_mot',0)
        self.outcrop_mot = kwargs.get('outcrop_mot',0)
        self.MRD=kwargs.get('MRD',np.zeros(len(t),))
#        self.AccTHo6 = kwargs.get('AccTHo6',np.zeros(len(t),))
#        self.AccTHoutcrop6 = kwargs.get('AccTHoutcrop6',np.ones(len(t),))
#        self.AccTHotype6 = kwargs.get('AccTHotype6',np.zeros(len(t),))
#        self.TauTH7a = kwargs.get('TauTH7a',np.zeros(len(t),))
#        self.GamTH7b = kwargs.get('GamTH7b',np.zeros(len(t),))
#        self.RespSpec9 = kwargs.get('RespSpec9',np.zeros(len(t),))
        self.OCR = kwargs.get('OCR',None)
        self.DS_soil = kwargs.get('DS_soil', None)
        self.DS_N = kwargs.get('DS_N', None)
        self.DS_freq = kwargs.get('DS_freq', None)
        self.filename = kwargs.get('filename')
        self.MRD = kwargs.get('MRD')
        self.D50= kwargs.get('D50')
        self.Cu= kwargs.get('Cu')
        self.gam_c= kwargs.get('gam_c')
        self.Su= kwargs.get('Su')
        self._calcInsituStress()


    def _calcInsituStress(self):
        '''
        Calculates the mid-depth of each layer, the vertical stress,
        the effective stresses, and the mean effective stresses.

        This method is automatically called by the soil profile
        object __init__.  For an external call, use SoilProfile.calcInsituStress.

        Notice that if ko is not given, it is calculated from phi.  If
        phi is not given, it is calculated from N160.  If none of ko, phi, or
        N160 are given, the ko value is assumed to be 0.5.  This is important
        because ko is used to calculate the mean effective stress which is
        used, in some cases, to calculate the degradation curves.

        Each of the inputs should be an array with the [0] value being the
        value of the layer at the surface:

            t = thickness of the layer

            unit_wt = the unit weight of the layer

            gwt = depth of the ground water table from the surface (positive value-
                    not an array)

            ko =  at rest earth pressure coefficient

            phi = effective friction angle.  Used to calc the ko if ko is given as
                    zero.

            N160 = the overburden- and energy-corrected SPT blowcount.  Used to
                    calc the ko it ko and phi are zero.
                    If N160, phi and ko are all zero for a layer, Ko is assumed to
                    be 0.5

        The following 'data attributes' or variables are created:
        self.t_mid
        self.sigv
        self.sigveff
        self.sigmeff
        '''
        t = self.t
        unit_wt = self.unit_wt
        gwt = self.gwt
        ko, phi, N160 = self.ko, self.phi, self.N160
        t_mid = np.zeros(len(t))
        sigv = np.zeros(len(t))
        sigveff = np.zeros(len(t))
        for i in range(len(t_mid)):
            if i == 0:
                t_mid[i] = (t[i] / 2)
                sigv[i] = (t_mid[i] * unit_wt[i])
            else:
                t_mid[i] = (t_mid[i - 1] + t[i - 1] / 2 + t[i] / 2)
                sigv[i] = (sigv[i - 1] + t[i - 1] / 2 * unit_wt[i - 1] +
                           t[i] / 2 * unit_wt[i])
            if t_mid[i] > gwt:
                sigveff[i] = (sigv[i] - 9.81 * (t_mid[i] - gwt))
            else:
                sigveff[i] = (sigv[i])
            if ko[i] == 0:
                if phi[i] != 0:
                    ko[i] = 1 - np.sin(phi[i] * np.pi / 180.)  # Using phi to calc ko
                elif N160[i] != 0:
                    phi[i] = (20 * N160[i]) ** 0.5 + 20  # Jaky, 1944, using N160 to calc phi and ko
                    ko[i] = 1 - np.sin(phi[i] * np.pi / 180.)
                    # Hatanae and Uchida, 1996
                else:
                    ko[i] = 0.5  # making an assumption
        sigmeff = sigveff * (1. + 2. * ko[:len(sigveff)]) / 3.
        self.t_mid, self.sigv = t_mid, sigv
        self.sigveff, self.sigmeff = sigveff, sigmeff


    def calcInsituStress(self, *args, **kwargs):
        '''
        This small function calculates the stresses at new depths
        Positional Arguments:
        Depths at which the insitu stresses should be calculated.
        If specified, the keyword argument "depths" is ignored.
        If neither the positional or keyword arguments are used,
        the mid-depth of each layer is used.

        Keyword Arguments:
        out = 'total':  the total vertical stress is output,
        Example:
            sigv = self.calcInsituStress([0.5,5.2], out='total')

        out = 'all':    the total, effective and mean effective stresses are output,
        Example:
             sigv, sigveff, sigmeff = self.calcInsituStress([0.5,5.2],
                                        out='all')

        out = 'eff':    the effective vertical stress is output.
        Example:
            sigveff = self.calcInsituStress([0.5,5.2], out='eff')

        out = 'meff':    the mean effective vertical stress is output.
        Example:
            sigmeff = self.calcInsituStress([0.5,5.2], out='meff')
        '''
        out = kwargs.get('out', 'all')
        gam_water =kwargs.get('gam_water', 9.81)
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)
        gwt = self.gwt
        unit_wt = self.unit_wt
        if np.isscalar(depths):
            depths = np.array([depths], dtype=float)
        try:
            depths = depths.ravel()
        except AttributeError:
            depths = np.array(depths, dtype=float)
            depths = depths.ravel()
        h = np.append(self.t,1000)
        hsum = np.cumsum(h, dtype=float) 
        sigv = np.zeros(len(depths))
        buoy = np.zeros(len(depths))
        ko = np.zeros(len(depths))

        # Find the layer for each depth and calc the total stresses
        j = 0
        for i in range(len(depths)):
            while depths[i] >= round(hsum[j],3):
                j += 1
            if depths[i] < h[0]:
                d = depths[i]
            else:
                d = depths[i] - hsum[j-1]

            if depths[i] >= gwt:
                buoy[i] = (depths[i] - gwt) * gam_water
            else:
                buoy[i] = 0

            if j == 0:
                sigv[i] = unit_wt[j] * d
            else:
                sigv[i] = np.sum(h[:j] * unit_wt[:j])+ unit_wt[j] * d
            ko[i] = self.ko[j]

        if out == 'total':
            return sigv
        elif out == 'all':
            sigveff = sigv - buoy
            sigmeff = sigveff * (1. + 2. * ko) / 3.
            return sigv, sigveff, sigmeff
        elif out == 'meff':
            sigveff = sigv - buoy
            sigmeff = sigveff * (1. + 2. * ko) / 3.
            return sigmeff
        else:
            sigveff = sigv - buoy
            return sigveff


    def calcVs30(self,**kwargs):
        '''
        Calculates the Vs30 of the profile according to eqn 20.4-1 of
        ASCE 7-10.


        Optional Keyword Argument:
            SiteClass = True or False [default].
                        If true, the site class according to ASCE 7-10 will
                        also be output with the Vs30.
                        If  False, the thickness of the entire soil profile
                        is returned instead.
            returns = True or False [default].  Whether or not anything is
                        returned, the following 'data attributes' or variables
                        are created:

                        self.Vs30
                        self.SiteClass

        Example:

            Vs30 = SoilObject.calcVs30(Vs_array, t_array, returns = True)
            Vs30,SiteClass = SoilObject.calcVs30(Vs_array,
                                                 t_array,
                                                 SiteClass=True,
                                                 returns = True)
        '''
        t = np.zeros(self.Vs.shape, dtype=float)
        t[:-1] = self.t
        Vsprofile = self.Vs
        returns = kwargs.get('returns',False)
        SiteClass = kwargs.get('SiteClass',False)

        n=-1
        if np.sum(t) <= 30:
            t[-1] = 30 - np.sum(t[:-1])
        else:

            t2 = np.zeros(t.shape,dtype=float)
            while np.sum(t2) < 30:
                n += 1
                t2[n] = t[n]

            if np.sum(t2) > 30:
                t2[n] = t2[n] - (np.sum(t2)-30)

            t = t2

        self.Vs30 = np.sum(t) / np.sum(t/Vsprofile)

        if SiteClass == True:
            if self.Vs30 > 1524:
                self.SiteClass = 'A'
            elif self.Vs30 >= 762:
                self.SiteClass = 'B'
            elif self.Vs30 >= 365.76:
                self.SiteClass = 'C'
            elif self.Vs30 >= 182.88:
                self.SiteClass = 'D'
            else:
                self.SiteClass = 'E'

            if returns == True:
                return self.Vs30, self.SiteClass
        else:
            if returns == True:
                return self.Vs30, np.sum(t)


    def calcAvgVs(self, depth):
        '''
        Calculates the average shear wave velocity of the top
        'depth' of the profile according to eqn 20.4-1 of
        ASCE 7-10.


        '''
        t = np.zeros(self.Vs.shape, dtype=float)
        t[:-1] = self.t
        Vsprofile = self.Vs
        n=-1
        if np.sum(t) <= depth:
            t[-1] = depth - np.sum(t[:-1])
        else:

            t2 = np.zeros(t.shape,dtype=float)
            while np.sum(t2) < depth:
                n += 1
                t2[n] = t[n]

            if np.sum(t2) > depth:
                t2[n] = t2[n] - (np.sum(t2)-depth)

            t = t2

        return np.sum(t) / np.sum(t/Vsprofile)


    def plot(self, *args, **kwargs):
        '''
        Outputs standard plots that illustrate the soil profile.

        Positional Arguments:

            List the plots you want to see.
            Possible plots:
                Vs
                unit_wt
                N160
                sigv
                sigmeff
                sigveff
                phi
                ko
                all

            Or, you can leave it blank and all will be output.

        Keywords Arguments:

            save = False (default) or True
                For now, if save = True, then the figures will also be closed.
                (Haven't figured it out yet.)

            extension = '.png' (default) or any other that matplotlib handles

            layerlabels = array of labels for each layer

            If layerlabels is specified, the following kwargs may be specified:
            labeldepths = depth at which the label should be placed
            labelxvals = x value at which the label is placed.

            Note that I have only implemented the labeling for the Vs plot.

        Example:

            self.plot('Vs','unit_wt', save=True, extension='.png')

            No data attributes are created with this method.
        '''
        if len(args) == 0:
            args = ('Vs','unit_wt','N160','sigv','sigmeff',
                        'sigveff','phi','ko')
        if 'all' in args:
            args = ('Vs','unit_wt','N160','sigv','sigmeff',
                        'sigveff','phi','ko')
        save = kwargs.get('save',False)
        ext = kwargs.get('extension','.png')
        lw = kwargs.get('lw',2)
        mark = kwargs.get('mark','-k')
        layerlabels = kwargs.get('layerlabels', False)
        if layerlabels:
            layername = kwargs.get('layername',self.layername)
            labeldepths = kwargs.get('labeldepths')
            labelxvals = kwargs.get('labelxvals')



        def profilecreator(array, **kwargs):
            '''
            Makes points at the top and bottom of each layer to make plotting
            better.
            '''
            depth = kwargs.get('depth',False)
            if depth == True:
                output = [0]
                for i in range(len(array)):
                    if output[i] == 0:
                        output.append(array[i])
                    else:
                        output.append(array[i-1])
                        output.append(array[i])
                output.append(array[-1])
                output.append(array[-1]*1.05)
            else:
                output = []
                for i in range(len(array)):
                    output.append(array[i])
                    output.append(array[i])
            return output

        def labeler(depths, labels, **kwargs):
            '''
            Adds labels to the plots.
            '''
            [xmin, xmax, ymin, ymax]  = plt.axis()
            xvals = kwargs.get('xvals',xmin)

            for i,z in enumerate(depths):
                plt.annotate(labels[i], xy=(xvals[i],z))


        depths = profilecreator(np.cumsum(self.t),depth=True)
        t_mid = self.t_mid

        if 'unit_wt' in args:
            unit_wt = profilecreator(self.unit_wt)
            plt.figure('unit_wt')
            plt.plot(unit_wt,depths,'-b')
            plt.xlabel('Unit Weight, $\gamma$  $(kN/m^3)$')
            plt.ylabel('Depth ($m$)')
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])

        if 'sigv' in args:
            sigv = self.sigv
            plt.figure('stress')
            plt.plot(sigv,t_mid,'-b',label='$\sigma_v$')
            plt.xlabel('Stress (kPa)')
            plt.ylabel('Depth ($m$)')
            plt.legend()
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])

        if 'sigmeff' in args:
            sigmeff = self.sigmeff
            plt.figure('stress')
            plt.plot(sigmeff,t_mid,'-r',label='$\sigma \'_m$')
            plt.xlabel('Stress (kPa)')
            plt.ylabel('Depth ($m$)')
            plt.legend()
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])

        if 'sigveff' in args:
            sigveff = self.sigveff
            plt.figure('stress')
            plt.plot(sigveff,t_mid,'-g',label='$\sigma \'_v$')
            plt.xlabel('Stress(kPa)')
            plt.ylabel('Depth ($m$)')
            plt.legend()
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])

        while 'phi' in args:
            phi = profilecreator(self.phi)
            if max(phi) == 0:
                break
            plt.figure('phi')
            plt.plot(phi,depths[:-2],'-b')
            plt.xlabel('Friction Angle, $\phi$')
            plt.ylabel('Depth ($m$)')
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            break

        if 'Vs' in args:
            Vs = profilecreator(self.Vs)
            plt.figure('Vs')
            plt.plot(Vs,depths,mark,lw=lw)
            plt.xlabel('Shear Wave Velocity, $V_s (m/s)$')
            plt.ylabel('Depth ($m$)')
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])

            if layerlabels:
                labeler(labeldepths, layername, xvals=labelxvals)


        while 'N160' in args:
            N160 = profilecreator(self.N160)
            if max(N160) == 0:
                break
            plt.figure('N160')
            plt.plot(N160,depths[:-2],'-b')
            plt.xlabel('Corrected SPT Blowcounts, $N_{1,60}$')
            plt.ylabel('Depth ($m$)')
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            break

        while 'ko' in args:
            ko = profilecreator(self.ko)
            if max(ko) == 0:
                break
            plt.figure('ko')
            plt.plot(ko,depths[:len(ko)],'-b')
            plt.xlabel('Lateral Earth Pressure Coefficient, $K_0$')
            plt.ylabel('Depth ($m$)')
            plt.title(self.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            break
        if save == True:
            import matplotlib
            figures=[manager.canvas.figure for manager in
                        matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
            for i, figure in enumerate(figures):
                figure.savefig(self.profile_name +'-'+figure._label+ext)
#                figure.clf()


    def calcDepth(self, stress, **kwargs):
        '''
        This method will return the depth at which a stress is felt.

        Keyword arguments:      Meaning/Options:
        stresstype                  'total', 'eff'[default], 'mean'
                                    corresponding to
                                    the total vertical stress, the vert
                                        effective stress, and the mean
                                        effective stress, respectively.
        Output:
        depths  = the depths at which the corresponding stresses are felt.
                    No data attributes are created with this method.

        Example:

            depths = self.calcDepth([50., 101.], stresstype='eff')
        '''
        unit_wt = self.unit_wt
        t = self.t
        t_mid = self.t_mid
        ko = self.ko
        gwt = self.gwt
        stresstype = kwargs.get('stresstype','eff')

        if np.isscalar(stress):
            stress = np.array([stress], dtype=float)
        try:
            stress.ravel()
        except AttributeError:
            stress = np.array(stress, dtype=float)
            stress.ravel()

        outdepths = np.zeros(stress.shape,dtype=float)

        if stresstype == 'total':
            stressprof = self.sigv
            unit_wt_eff = unit_wt[:-1]
        elif stresstype == 'eff':
            stressprof = self.sigveff
            unit_wt_eff = np.where(np.cumsum(t) <= gwt, unit_wt[:-1],
                                   unit_wt[:-1]-9.81)
        elif stresstype == 'mean':
            stressprof = self.sigmeff
            unit_wt_eff = np.where(np.cumsum(t) <= gwt, unit_wt[:-1],
                                   unit_wt[:-1]-9.81) * (1 + 2 * ko[:-1]) / 3
        else:
            print('I didn\'t understand the stresstype')
            return np.ones(stress.shape,dtype=float) * -1000.

        for i,sig in enumerate(stress):
            try:
                no = np.where(sig < stressprof)[0][0]
            except IndexError:
                print('Stress of {} kPa not reached in this profile.'.format(
                        sig))
                outdepths[i] = -1000
                continue


            stressinterface = stressprof[no-1] + (unit_wt_eff[no-1] *
                            t[no-1] / 2)
            if stressinterface >= sig:
                outdepths[i] = t_mid[no-1] + ((sig - stressprof[no-1]) /
                            unit_wt_eff[no-1])
            else:
                outdepths[i] = np.sum(t[:no]) + ((sig -
                                stressinterface) / unit_wt_eff[no])

        return outdepths




    def getLayer(self, depth):
        '''
        Return the layer number for a given depth.
        Example:

            layernums = self.getLayer(10.25)

            layernums = self.getLayer([10.25, 13.7])

            No data attributes are created with this method.
        '''
        tsum = np.cumsum(self.t)
        if np.isscalar(depth):
            depth = np.array([depth], dtype=float)
        try:
            depth.ravel()
        except AttributeError:
            depth = np.array(depth, dtype=float)
            depth.ravel()
        layernumbers = np.ones(depth.shape, dtype=int) * -1

        for i,z in enumerate(depth):
            layernumbers[i] = np.where(z <= tsum)[0][0]

        if len(layernumbers) == 1:
            return layernumbers[0]
        else:
            return layernumbers





class GroundMotion(object):
    '''
    This object acts as a container for all ground motion related variables.

    Keyword Arguments:          Definition/defaults:
        name           The name of the ground motion/ground motion file.

        filename       The filename of ground motion file.  If not specified,
                           name needs to be a valid filename in the path.

        scalefactor    Factor by which the motion time history is scaled.
                            This scaling occurs when the motion is read from
                            the file. If the motion is specified as a keyword
                            argument, no scaling occurs.

        Butrwrth       If True, a Butterworth filter is applied when the motion
                        is called from the file.  If the motion is specified
                        as a keyword argument, the Butterworth filter is not
                        applied. Default is False

        b_order        The order of the Butterworth filter.  Default is 4.

        cutofffreq     Cutoff frequency (in Hz) for the Butterworth filter.
                        Default is 100 Hz.

        g              Acceleration due to gravity. Default is 9.81.  This is
                        used to convert the motion from units of g to m/s/s.
                        Apply the appropiate value for the units of the motion
                         and the soil profile. This code assumes everything
                         is in SI units unless specified.

    ====  The following are only needed if a valid filename is not given. ====
       motion          The input acceleration time history. If this is not
                           specified, a valid filename must be given.

       NFFT            The number of points used in the FFT.  Only required if
                           the motion keyword argument is specified.

       FA              The Fourier Amplitude spectrum of the input motion.
                           Only required if the motion keyword argument is
                           specified.

       f                An array of frequencies corresponding to the FA in Hz.
                           Only required if the motion keyword argument is
                           specified.

       w               An array of angular frequencies corresponding to the FA
                           in rad/sec. Only required if the motion keyword
                           argument is specified.

       NPTS            Number of points in the input time history. Only
                           required if the motion keyword argument is specified.

       dt              The sample spacing of the time history input motion in
                           seconds. Only required if the motion keyword
                           argument is specified.




       dummy            Let dummy=True if you just want to calc the transfer function of the profile.


       'Data attributes' or variables created by the __init__ method:
       self.name
       self.filename
       self.scalefactor
       self.Butrwrth
       self.b_order (if self.Butrwrth)
       self.cutofffreq (if self.Butrwrth)
       self.g
       self.motion
       self.NFFT
       self.FA
       self.f
       self.w
       self.NPTS
       self.dt
       self.dummy
       self.time

    '''
    def __init__(self,**kwargs):
        self.name = kwargs.get('name')
        self.filename = kwargs.get('filename')
        formatin = kwargs.get('formatin','PEER')
        self.scalefactor = kwargs.get('scalefactor',1.)
        self.Butrwrth = kwargs.get('Butrwrth',False)
        if self.Butrwrth:
            self.b_order = kwargs.get('b_order',4.)
            self.cutofffreq = kwargs.get('cutofffreq',100)
        self.g = kwargs.get('g',9.81)
        self.motion = kwargs.get('motion')
        self.NFFT = kwargs.get('NFFT')
        self.FA = kwargs.get('FA')
        self.f = kwargs.get('f')
        self.w = kwargs.get('w')
        self.NPTS = kwargs.get('NPTS')
        self.dt = kwargs.get('dt')
        self.dummy = kwargs.get('dummy',False)
        if self.NPTS and self.dt:
            self.time = np.linspace(0,(self.NPTS-1) * self.dt,self.NPTS)
        if self.dummy == True:
            self.name = 'dummy'
            self.motion = np.ones([8192,],dtype=float) * -1.
            self.NPTS = 8192
            if self.NFFT == None:
                self.NFFT = 8192
            if self.dt == None:
                self.dt = 0.01
            self.f = np.linspace(0, 1/self.dt/2, self.NFFT//2+1)
            self.w = self.f[:self.NFFT//2+1] * 2 * np.pi
            self.FA = np.ones(self.w.shape, dtype='complex128')
        if self.motion == None:
            if self.filename == None:
                self.getInMotion(self.name, **kwargs)
            else:
                self.getInMotion(self.filename, **kwargs)

        if self.FA is None:
            self.getMotInfo(self.motion, self.dt, **kwargs)



    def getInMotion(self, motion_name,  **kwargs):
        '''
        Reads in the input motion from a PEER-formatted file.  Also converts it
        to the frequency domain and applies scaling.  Also, applies a Butterworth
        lowpass filter, if specified.

        Keyword Arguments:

        formatin = 'PEER' [default], 'PEER-OLD', or 'PEER-scaled'

        directout = False [default] or True. If true, the __init__ data
                    attributes are not updated.

        The following 'data attributes' or variables may be created:
        self.Ar  =  a normalized rms acceleration
        self.predfreq =  an estimation of the dominant frequency of the motion
        self.Ar2       self.Ar / self.predfreq
        '''
        Butrwrth = kwargs.get('Butrwrth', False)
        if self.Butrwrth:
            self.b_order = kwargs.get('b_order',4.)
            self.cutofffreq = kwargs.get('cutofffreq',100)
        g = self.g
        formatin = kwargs.get('formatin', 'PEER')
#        baseline = kwargs.get('baseline', False)
        extras = kwargs.get('extras', True)
        directout = kwargs.get('directout', False)
        # Import the ground motion
        count = 0
        eqmot = []
        f = open(motion_name, 'rb')
        if formatin == 'PEER-scaled':
           for line in f:
                count += 1
                if count == 7:
                    NPTS = int(line.split()[1][:-1])
                    dt = float(line.split()[3][:-4])
                if count > 7:
                    for num in line.split():
                        eqmot.append(float(num))
        elif formatin == 'PEER-OLD':
            for line in f:
                count += 1
                if count == 4:
                    line=str(line)
                    NPTS = int(line.split()[1].replace(',',''))
                    dt = float(line.split()[3])
                if count > 4:
                    for num in line.split():
                        eqmot.append(float(num))
            #Added in by Tat as part of Tyler's request
        elif formatin.lower() == "custom":      #This gm file is basically produced using Excel
            NPTS = 0
            time = []
            for line in f:
                if NPTS == 0:
                    time.append(float(line.split()[0]))
                if NPTS == 1:
                    time.append(float(line.split()[0]))   
                eqmot.append(float(line.split()[1]))
                NPTS += 1
            dt = time[1]-time[0]
        else:
            for line in f:
                count += 1
                if count == 4:
                    NPTS = int(line.split()[0])
                    dt = float(line.split()[1])
                if count > 4:
                    for num in line.split():
                        eqmot.append(float(num))
        motion = np.array(eqmot) * self.scalefactor
        if max(motion.shape) != NPTS :
            print('NPTS discrepancy! {}  {}'.format(max(motion.shape), NPTS))
            print(motion_name)
            NPTS = max(motion.shape)

        motion = motion.reshape([NPTS])
#        if baseline == True:
#            motion = motion - np.average(motion)
        f.close()
        del eqmot, count
        try:
            del num
        except UnboundLocalError:
            del time

        self.getMotInfo(motion, dt, **kwargs)


    def getMotInfo(self, scaledmotion, dt, **kwargs):
        '''
        This will get the transform and other info of a motion if only the
        acceleration time history is input.
        '''
        extras = kwargs.get('extras',True)
        Butrwrth = kwargs.get('Butrwrth', False)
        if self.Butrwrth:
            self.b_order = kwargs.get('b_order',4.)
            self.cutofffreq = kwargs.get('cutofffreq',100)
        g = self.g
        directout = kwargs.get('directout', False)
        g = self.g
        NPTS = len(scaledmotion)
        NFFT = 2 ** nextpow2(NPTS)
        f = np.linspace(0, 1/dt/2, NFFT//2+1)
        w = f[:NFFT//2+1] * 2 * np.pi
        if Butrwrth == True:
            w_cutoff = self.cutofffreq * 2 * np.pi
            butter = (1 / (1 + (np.abs(w) / w_cutoff) ** (2 * self.b_order))
                                                        ) ** 0.5
        else:
            butter =  np.ones(w.shape)

        FA = np.fft.rfft(scaledmotion  * g, NFFT) * butter


        if directout == True:
            return scaledmotion, NFFT, FA, f, w, NPTS, dt
            pass
        else:
            self.motion = np.fft.irfft(FA / g)[:NPTS]
            self.NFFT = NFFT
            self.FA = FA
            self.f = f
            self.w = w
            self.NPTS = NPTS
            self.dt = dt
            self.time = np.linspace(0,(NPTS-1) * dt,NPTS)
            if extras:
                sortedfreq = np.array(sorted(zip(np.abs(FA),f)))
                predfreq = (np.sum(sortedfreq[-NFFT//200:,0] *
                    sortedfreq[-NFFT//200:,1]) /
                    np.sum(sortedfreq[-NFFT//200:,0]))
                Ar = calcAriasRatio(scaledmotion, dt)
                Ar2 = Ar / predfreq
                self.Ar = Ar
                self.Ar2 = Ar2
                self.predfreq = predfreq




class Dynamics(object):
    '''
    This class combines the soil profile and the ground motion objects to
    produce the dynamic Equivalent Linear Procedure.

    Input:
        SoilObject   =   This is the object created by instantiating the
                            SoilProfile class

        MotionObject =   This is the object created by instantiating the
                            GroundMotion class.

            You can specify these two objects so that you can perform different
            analyses on the same profile with the same motion


        Keyword Arguments:          Definition and default values:
            modtype     =    specifies complex shear modulus calculations.
                                The options are:
                                'FreqInd': Frequency Independent Complex Shear
                                            Modulus (Kramer, 1996) [default]
                                'FreqDep': Frequency Dependent Complex Shear
                                            Modulus (Udaka, 1975)
                                            (used by SHAKE?)
                                'Simplified': Simplified complex shear modulus
                                (Kramer, 1996)
                The calculations can be found under the _calcWaves function.

            InitGratio       The initial Gratio estimate.  Should be either a
                                scalar or the length of the layers (not
                                counting bedrock). By default it is equal to
                                Gmax.

            BaseGratio       The initial Gratio for the bedrock.  Since
                                there are no iteration or degradation
                                curves for the bedrock, this should be 1. [Default]

            Dminfactor      Correction of Dmin. Default is one
                            

            BaseDamping     The damping constant of the bedrock, as a decimal.
                                0.01 is default. I may try to add the option of
                                adding degradation curves later

            DegradCurves    Specifies the method used to calculate the shear
                                modulus and damping degradation curves.
                                The only options currently supported are
                                    'IZ' Ishibashi and Zhang (1993) and
                                    'DS' for Darendeli and Stokoe (2001)
                                It shouldn't be difficult to add more.
                                If the 'DS' option is used, the DSVariables
                                object is instantiated and used.

           iters            Max number of iterations, default is 20

           Error            The threshold error (decimal, as a ratio of Gmax)
                               at which iterations will stop, default is 0.02.
           MRD_rand         Randomizes modulus reduction and damping curves
           
           epsilon          If you input this vector the epsilon values for
                            MRD curves are calculated by this amount deviation

  ===========>  Iterations continue until either the number of iterations is
                    reached or the maximum computed error of the iteration is
                    less than 'Error,' the threshold error.

           outcrop_mot      0 if the input motion is an outcrop, 1 otherwise.
                               This overrides the SoilProfile object data
                               attribute, if specified, otherwise, the
                               SoilObject value is used.
           loc_mot           Overrides the location of motion given in
                               the SoilProfile object data
                               attribute, if specified, otherwise, the
                               SoilObject value is used.

           strainratio      Ratio of equivalent uniform strain divided by the
                               maximum strain. This should be a single value;
                               the default is given by the SoilObject.
                               The default given by the SoilObject is 0.65.
                                If anything other than a numerical value is
                                given, a strain ratio is calculated using
                                dissipated energy to give a weighted average
                                of the strain time history.
                                This is experimental; see self.calcDissEn to
                                see how it works.

           start_iters      If False, code will wait until the
                                   DynamicsObject.iterate() method is called.
                                   Default is True.

           run_all          If True, the number of equivalent cycles and
                                dissipated energy will be calculated for all
                                mid-layer depths.

          verbose           True or False [default].  If True, details from
                                  each iteration are given.

========  The following options allow the user to finely specify degradation
            curves. They are all completely optional. They are useful when
            comparing outputs from different codes.   ========================

           DCO_sigmeff      Overrides the calculated mean effective stresses
                               for use in the degradation curve equations.
                               Default are the SoilObject values.

           DCO_PI           Overrides soil profile PI values for use
                               in the degradation curve equations.
                               Default are the SoilObject values.

           DS_OCR           Provides OCR values for the Darendeli and
                               Stokoe (2001) Degradation Curves.
                               Default is 1 for all layers.

           DS_soil          Soil number for Darendeli and Stokoe Degrad Curves
                               The soil type indicator:
                                   0 for a general set of curves
                                     (Table 8.12 from Darendeli's Dissertation)
                            		 1 for clean sand
                                   2 for sands with high fines contents
                                   3 for silts
                         		 4 for clays

           DS_N             Number of cycles for the Darendeli and
                               Stokoe (2001) Degradation Curves.
                               Default is 10 for all layers.

           DS_freq          Frequency for the Darendeli and
                               Stokoe (2001) Degradation Curves.
                               Default is 1 Hz for all layers.

           DS_Dmin           An alternative way to calc Dmin for the Darendeli
                               and Stokoe Degradation curves.
                               1 for 'Default', 0 for 'Green'
                               !!! Not Implemented !!!

           DCO_curve        If True, user must input their own degradation
                               curves, specified with the following keywords:
                DCO_gam = shear strains as decimals (an array n long)
                DCO_Gratio = G/Gmax values (decimal) (an array m  x n, where
                            m is the number of layers above the baserock)
                DCO_damping = damping (decimal) in the same shape as DCO_Gratio

       The following data attributes are created by the __init__ method.
           Since the init method calls other methods, the list may not be
           complete.

           self.g
           self.Gmax
           self.G
           self.D
           self.DegradCurves
           self.modtype
           self.iters
           self.ComIter     =  Flag to show when iterations have completed.
           self.Error
           self.gam
           self.Gratio
           self.damping
           self.loc_mot
           self.SoilObj
           self.MotObj
           self.outcrop_mot
           self.dissEn
           self.strainratio
           self.verbose
           self.Ar  =  a normalized rms acceleration.
           self.Ar_ratio = self.Ar / Ar of the input motion.


    '''
    def __init__(self, SoilObject, MotionObject, **kwargs):
        self.g = SoilObject.g
        self.Gmax = (SoilObject.unit_wt / self.g) * SoilObject.Vs ** 2
        self.Dmin = np.copy(self.Gmax[:-1])
        self.G = np.copy(self.Gmax)
        self.G[:-1] = self.G[:-1] * kwargs.get('InitGratio', 1.)
        self.G[-1] = self.G[-1] * kwargs.get('BaseGratio',1.)
        self.ComIter = False
        self.D = np.ones(len(self.Gmax)) * 0.05
        self.D[-1] = kwargs.get('BaseDamping',0.01)  # <<<---------------------------
        self.modtype = kwargs.get('modtype', 'FreqInd')
        self.iters = kwargs.get('iters',30)
        self.Error = kwargs.get('Error',0.02)
        # Version 1
        self.Dminfactor = kwargs.get('Dminfactor', 1)
        self.MRD_rand = kwargs.get('MRD_rand',0)
        self.epsilon = kwargs.get('epsilon',np.ones(len(self.Gmax[:-1])) * -10)
        gam = np.logspace(-7, -1.2, num=25, endpoint=True)
        gam_n = np.ones(len(gam))
        self.Gratio = np.zeros([len(SoilObject.PI),len(gam)],dtype=float)
        self.damping = np.zeros([len(SoilObject.PI),len(gam)],dtype=float)
        # Get the degradation curves
        if kwargs.get('DCO_curve',False) == True:
                self.gam = kwargs.get('DS_gam')
                self.Gratio = kwargs.get('DCO_Gratio')
                self.damping = kwargs.get('DCO_damping')
        for count in range(len(SoilObject.PI)):
            gam_c = SoilObject.gam_c[count] # The cut off strain
            for i in range(len(gam)):
                if gam[i]>gam_c :
                   gam_n[i]=gam_c
                else:
                   gam_n[i]=gam[i]
            if SoilObject.MRD[count] == 'DS':
                sigmeff = kwargs.get('DCO_sigmeff',SoilObject.sigmeff[count])
                PI= kwargs.get('DCO_PI',SoilObject.PI[count])
                if len(SoilObject.OCR) == 0:
                    OCR= kwargs.get('DS_OCR',1) #Non-over consolidated as default
                else:
                    OCR = SoilObject.OCR[count]
                soiltype = kwargs.get('DS_soil', SoilObject.DS_soil[count])
                if soiltype == None:
                    soiltype = 0 #General default
                N = kwargs.get('DS_N', SoilObject.DS_N[count])
                if N == None:
                    N = 10 # 10 cycles default
                freq = kwargs.get('DS_freq',SoilObject.DS_freq[count])
                if freq == None:
                    freq = 1  # 1 Hz default
    #           self.DSvars = DSvariables(kwargs,SoilObject)
                self.gam,self.Gratio[count], self.damping[count], self.Dmin[count] = dc.DS_2001(gam=gam_n,
                                                      PI=PI,
                                                      sigm=sigmeff,
                                                      OCR=OCR,
                                                      soil=soiltype,
                                                      N=N,
                                                      frq=freq,
                                                      MRD_rand=self.MRD_rand,
                                                      epsilon=self.epsilon[count],
                                                      Su=SoilObject.Su[count],
                                                      Gmax=self.Gmax[count])

            elif SoilObject.MRD[count] == 'IZ':
                self.gam,self.Gratio[count], self.damping[count], self.Dmin[count] = dc.IZ_1993(gam = gam_n,
                                                        PI=SoilObject.PI[count],
                                                        sigm=SoilObject.sigmeff[count])

            elif SoilObject.MRD[count] == 'Menq':
                self.gam,self.Gratio[count], self.damping[count], self.Dmin[count] = dc.Menq(gam = gam_n,
                                                        Cu = SoilObject.Cu[count],
                                                        D50 = SoilObject.D50[count],
                                                        frq = SoilObject.DS_freq[count],
                                                        N = SoilObject.DS_N[count],
                                                        sigm = SoilObject.sigmeff[count],
                                                        MRD_rand = self.MRD_rand,
                                                        epsilon = self.epsilon[count],
                                                        Su = SoilObject.Su[count],
                                                        Gmax = self.Gmax[count])

            elif SoilObject.MRD[count] == 'Vucetic_Dobry':
                self.gam,self.Gratio[count], self.damping[count] = dc.Vucetic_Dobry (gam = gam_n, PI=SoilObject.PI[count])
                self.Dmin[count]=0.01
            elif SoilObject.MRD[count] == 'PietPeat':
                self.gam,self.Gratio[count], self.damping[count], self.Dmin[count] = dc.PietPeat(gam = gam_n,
                                                                               sigm = SoilObject.sigmeff[count],
                                                                               N = SoilObject.DS_N[count],
                                                                               frq = SoilObject.DS_freq[count],
                                                                               MRD_rand = self.MRD_rand,
                                                                               epsilon = self.epsilon[count],
                                                                               Su = SoilObject.Su[count],
                                                                               Gmax = self.Gmax[count])
            elif SoilObject.MRD[count] == 'Linear':
                self.gam = gam
                self.Gratio[count] =  np.ones(len(self.gam))
                self.damping[count] = np.ones(len(self.gam)) * 0.005/self.Dminfactor # The final value of rock damping should be 0.005
                self.Dmin[count]=0.005/self.Dminfactor
            elif SoilObject.MRD[count] == 'Undamped_linear':
                self.gam = gam
                self.Gratio[count] =  np.ones(len(self.gam))
                self.damping[count] = np.zeros(len(self.gam))
            
            else:
                line=('Degradation Method input not understood\n' +
                        'or a custom degradation curve has been entered.')
            
            self.damping[count]=self.damping[count] + (self.Dminfactor-1) * self.Dmin[count]
            self.gam=gam
            
        self.D[:-1] = [damping_layer[0] for damping_layer in self.damping]
        #End
        self.SoilObj = SoilObject
        self.loc_mot = kwargs.get('loc_mot', self.SoilObj.loc_mot)

        self.MotObj = MotionObject
        self.outcrop_mot = kwargs.get('outcrop_mot',self.SoilObj.outcrop_mot)
        self.gammf = None
        self.dissEn = []   #previously saved as None, and later compared between an [] and a None. Works for Python2 apparently

        # Start Iterations
        start_iters = kwargs.get('start_iters',True)
        if start_iters:
            self.strainratio = kwargs.get('strainratio',self.SoilObj.strainratio)
            self.verbose = kwargs.get('verbose',False)
            self.iterate(iters=self.iters,
                         Error=self.Error,
                         strainratio=self.strainratio,
                         verbose=self.verbose)

        # Run_all
        run_all = kwargs.get('run_all',False)
        if run_all:
            self.calcNeq()
            self.Ar,self.PeakTimeRatio = calcAriasRatio(self.calcAcc(domain='time',
                                          MotType='within',
                                          MaxOnly=False,
                                          ), self.MotObj.dt, returnall=True)
            Ar_input = calcAriasRatio(self.MotObj.motion,self.MotObj.dt)
            self.Ar_ratio = self.Ar / Ar_input


    def iterate(self,**kwargs):
        '''
        Performs the equivalent linear iterations until either the iteration
        limit is reached or the error goes below the predefined threshold.

        Keywords:             Definition/Default:
          iters                number of iterations, default is that defined
                                  by the self.__init__ function

          Error                the error threshold, default is that defined
                                  by the self.__init__ function

          strainratio          Ratio of equivalent uniform strain divided by
                                  the maximum strain. This should be a single
                                  value; the default is given
                                  by the SoilObject.  The default given by the
                                  SoilObject is 0.65.
                                  If anything other than a numerical value is
                                  given, a strain ratio is calculated using
                                  dissipated energy to give a weighted average
                                  of the strain time history.

          verbose               True or False [default].  If True, details from
                                  each iteration are given.

        Data attributes created by this method:
        self.gamavgarray = average shear strain of all layers
        self.Garray = array of current degraded shear moduli for every layer
                        for every iteration
        self.Darray = array of current degraded damping for every layer
                        for every iteration
        self.Errarray = array of errors for every layer
                        for every iteration
        self.Gratiocurv = an object to interpolate new values of modulus
        self.Dcurv = an object to interpolate new values of damping
        self.count = current or final count of iterations
        self.ratio = strain ratios used for the iteration
        self.ReportError = 'y' if an error has occurred, 'n' if not.

        Updated or created data attributes:
        self.G
        self.D
        self.A
        self.B
        self.ComIter
        self.gammf = strain for mid point or each layer in freq domain
        self.taumf = stress for mid point or each layer in freq domain
        '''
        self.gamavgarray = np.zeros((self.Gratio.shape[0],max(self.iters,1)),
                                    dtype=float)
        self.Garray = np.zeros((self.Gratio.shape[0],max(self.iters,1)),
                               dtype=float) #G/Gmax
        self.Darray = np.zeros((self.Gratio.shape[0],max(self.iters,1)),
                               dtype=float) #D
        self.Errarray = np.zeros((max(self.iters,1), len(self.G)-1),
                dtype=float)
#        self.Gratiocurv = interp1d(self.gam, self.Gratio, kind='linear', axis = 1)
#        self.Dcurv = interp1d(self.gam, self.damping, kind='linear', axis = 1)
        self.Gratiocurv = interp1d(np.log10(self.gam), self.Gratio,
                                   kind='linear', axis = 1)
        self.Dcurv = interp1d(np.log10(self.gam), self.damping,
                              kind='linear', axis = 1)
        self.iters = kwargs.get('iters', self.iters)
        Error = kwargs.get('Error', self.Error)
        self.strainratio = kwargs.get('strainratio',self.SoilObj.strainratio)
        self.verbose = kwargs.get('verbose',False)
        iterate = True
        self.count = 0
        while iterate == True:  # Begin Iterations
            self.count += 1
            self._calcWaves()  # Get the new upgoing and down going waves
            try:
                # Won't work if the strain ratio is a string, go to except
                self.ratio = (float(self.strainratio) *
                                np.ones(self.G[:-1].shape))
                # Get the max stress and strain values
                taumax,gammax = self._calcTauGam(MaxOnly=True,domain='time',
                                             returns=True)

            except ValueError:
                # Get the stress and strain for each layer in the time domain
                taut, gamt = self._calcTauGam(MaxOnly=False,domain='time',
                                          returns=True)
                # Calc the strain ratio using the dissipated energy method
                self.ratio = self.calcStrainRatio_DE(taut, gamt)
                gammax = np.max(np.abs(gamt),axis=1)

            self.gamavgarray[:,self.count-1] = self.ratio * gammax
            self.gammax = gammax
            # For really stiff profiles or soft motions, self.ratio * gammax
            # may fall below the gam range. Here we make sure we can interpolate.
            gammain = np.min((
                            np.max((
                                self.ratio*gammax,
                                self.gam[0] * np.ones(
                                gammax.shape)), axis=0),
                                self.gam[-1]*np.ones(gammax.shape)
                              ), axis=0
                            )
            # Get the new G and D from interpolated curves
            try:
                Gnewmat = self.Gratiocurv(np.log10(gammain))
#                Gnewmat = self.Gratiocurv(self.ratio * gammax)
                Gnew = self.Gmax[:-1] * np.diag(Gnewmat)
                Dnewmat = self.Dcurv(np.log10(gammain))
#                Dnewmat = self.Dcurv(self.ratio * gammax)
                Dnew = np.diag(Dnewmat)
            except ValueError:
                print('ValueError!')
                print('Soil Profile: {}'.format(self.SoilObj.profile_name))
                print('Ground Motion:  {}  {}'.format(self.MotObj.name,
                            self.MotObj.filename))
                print('Strain Ratio: {}'.format(self.ratio))
                print('Gamma Max: {}'.format(gammax))
                print('G: {}'.format(self.G))
                print('D: {}'.format(self.D))
                self.ReportError = 'y'
                pdb.set_trace()
                raise Exception("ValueError in the interpolation of the" +
                               'degradation curves'  )
                return

            Err = (self.G[:-1] - Gnew)/self.Gmax[:-1]  # calc the error
            self.Errarray[self.count - 1,:] = Err
            # Check if we can duck out of iterations
            if np.max(abs(Err)) <= Error:
                iterate = False
                self.ReportError = 'n'
                self.ComIter = True
                self.MaxError = np.max(abs(Err))
            if self.count >= self.iters:
                iterate = False
                if np.max(abs(Err)) <= Error:
                    pass
                    self.ReportError = 'n'
                    self.ComIter = True
                else:
                    print('Reached {:d} iterations without converging!'.format(
                                            self.iters))
                    self.ReportError = 'y'
                    self.ComIter = False
                self.MaxError = np.max(abs(Err))
            # print out data about each iteration/layer
            if self.verbose == True:
                print('Iteration: {}\n'.format(self.count))
#                print('Gamma Max: {}'.format(gammax))
#                print('G:'.format(self.G))
#                print('D:'.format(self.D))
#                print(np.max(self.A))
                print('{:^11} {:^9} {:^9} {:^9} {:^9} {:^9}'.format(
                            'Gamma_max:',
                            'Gold',
                            'Gnew',
                            'Dold',
                            'Dnew',
                            'Error'))
                for i in range(len(Gnew)):
                    print('{:^10.6e} {:^9.0f} {:^9.0f}'.format(
                            gammax[i],
                            self.G[i],
                            Gnew[i]) +
                          '{:^9.4f} {:^9.4f} {:^9.4}'.format(
                            self.D[i],
                            Dnew[i],
                            Err[i]))
            # Get ready for the new iteration.
            if iterate == True:
                self.G[:-1] = Gnew
                self.D[:-1] = Dnew
                self.Garray[:,self.count-1] =  self.G[:-1]/self.Gmax[:-1]
                self.Darray[:,self.count-1] = self.D[:-1]



    def _calcWaves(self):
        '''
        Calculates the up-going and down-going wave amplitudes.

            This method is called without input arguments to calculate the
            amplitudes of the upgoing and downgoing waves (A & B).  The input
            are already defined as part of the object.

            The following data attributes are created:
            self.Gstr = complex shear modulus
            self.Vsstr = complex shear wave velocity
            self.ks = wave number?
            self.alps = alpha, the layer impedence ratio
            self.A = amplitude of upgoing waves
            self.B = amplitude of downgoing waves
            self.Updated = flag for something.
        '''
        g = self.g
        G = self.G
        D = self.D
        h = self.SoilObj.t
        unit_wt = self.SoilObj.unit_wt
        try:
            w = self.MotObj.w
            FA = self.MotObj.FA
        except AttributeError:
            print(self.MotObj.filename)
            print(self.MotObj.dt)
        # From the Deepsoil Manual:
        if self.modtype == 'FreqDep':  # Frequency Dependent Complex Shear Modulus (Udaka, 1975) (Shake?)
            self.Gstr = G * (1 - 2 * D ** 2 + 1j * 2 * D * (1 - D ** 2) ** (0.5 + 0j))
        elif self.modtype == 'FreqInd':  # Frequency Independent Complex Shear Modulus (Kramer, 1996)
            self.Gstr = G * (1 + 1j * 2 * D)
        elif self.modtype == 'Simplified':  # Simplified complex shear modulus (Kramer, 1996)
            self.Gstr = G * (1 - D ** 2 + 1j * 2 * D)
        Gstr = self.Gstr
        self.Vsstr = (Gstr * g / unit_wt) ** 0.5
#        self.Vsstr = ne.evaluate('(Gstr * g / unit_wt) ** 0.5')
        Vsstr = self.Vsstr
        self.ks = w / Vsstr[:,np.newaxis]
        ks = self.ks

        alps = ((unit_wt[:-1] * Vsstr[:-1]) /
                (unit_wt[1:] * Vsstr[1:]))
        # Initialize the up and down going waves.
        A=np.ones(ks.shape, dtype='complex128')
        B=np.ones(ks.shape, dtype='complex128')
        # Calculate the waves, this is somewhat computationally intensive.  C code could be useful here in the future.
        for i in range(len(h)):
            A[i+1,:] = (0.5 * A[i,:] * (1 + alps[i]) *
                    np.exp(1j * ks[i,:] * h[i]) +
                    0.5 * B[i,:] * (1 - alps[i]) *
                    np.exp(-1j * ks[i,:] * h[i]))
            B[i+1,:] = (0.5 * A[i,:] * (1 - alps[i]) *
                    np.exp(1j * ks[i,:] * h[i]) +
                    0.5 * B[i,:] * (1 + alps[i]) *
                    np.exp(-1j * ks[i,:] * h[i]))
        # Get the input motion
        FD = np.zeros(FA.shape, dtype=complex)
        FD[1:] = FA[1:] / (1j*w[1:]) ** 2
        FD[0] = FA[0] #* (self.MotObj.dt * (self.MotObj.NPTS - 1)) ** 2   # <<<<==================== This may not be strictly legitimate
        if self.loc_mot == 0:
            index = -1
        else:
            index = self.loc_mot - 1   # Subtracting one here to convert to python's 0 counting method
        if index > len(A):
            index = -1
        if self.outcrop_mot == 0:
            A1new = FD / (2 * A[index, :])
        else:
            A1new = FD / (A[index,:] + B[index, :])

#        Aratio = A1new / A[index,:]    # <++++  Double checking my calcs: 3 Sept 2013, Everything is okay.
#        A2 = A * Aratio[np.newaxis,:]
#        Bratio = A2[0,:] / A[0,:]
#        B2 = B * Bratio[np.newaxis,:]
        # Getting the corrected waves from the input motion
        self.A = A * A1new[np.newaxis, :]
        self.B = B * A1new[np.newaxis, :]
#        plt.figure(1)
#        plt.plot(abs(self.A[0,:]))
#        plt.plot(abs(A2[0,:]),'o')
#        plt.figure(2)
#        plt.plot(abs(self.B[0,:]))
#        plt.plot(abs(B2[0,:]),'o')
#        plt.show()
#        pdb.set_trace()
        self.alps = alps
        self.Updated = True


    def _getDegrad(self,Degrad, PI, sigmeff, **kwargs):
        '''
        Gets the shear modulus and damping degradation curves.

        Outputs a 1D array (gam = shear strain (decimal) and 2- 2D arrays
        (Gratio = G/Gmax and damping (decimal)) that represent the shear
        modulus and damping degradation curves.  It depends on the module
        degrcurv, so use that if you want to customize it.  Or, create your own
        degradation curves.
        Right now, Degrad must equal 'IZ' for the Ishibashi and Zhang (1993)
        curves  or 'DS' for Darendeli and Stokoe (2001).
        sigmeff must be in kPa

        Keyword Arguments:
            (See the __init__ for these definitions, only for 'DS')
        OCR
        soil
        N
        freq
        Dmintype !!! Not Implemented !!!

        Everything but Degrad should be an array, or you will get errors.
        '''
        #straindef = np.array([0.0000001, 0.000001, 0.000003, 0.00001,
                                #0.00003, 0.0001, 0.0003, 0.001, 0.003,
                                #0.01, 0.99])
        straindef = np.logspace(-6, -2, num=20, endpoint=True)
        if Degrad == 'DS':
            OCR = kwargs.get('OCR')
            soil = kwargs.get('soil')
            N = kwargs.get('N')
            freq = kwargs.get('freq')
#            Dmintype = kwargs.get('DS_Dmin')
            if np.isscalar(OCR):
                OCR = np.ones(len(PI),dtype=float) * OCR
            if np.isscalar(soil):
                soil = np.ones(len(PI),dtype=int) * soil
            if np.isscalar(N):
                N = np.ones(len(PI),dtype=float) * N
            if np.isscalar(freq):
                freq = np.ones(len(PI),dtype=float) * freq
        gam = kwargs.get('gam', straindef)
        if gam == None:
            gam = straindef
        damping = np.zeros([len(PI),len(gam)],dtype=float)
        Gratio = np.zeros([len(PI),len(gam)],dtype=float)
        Dmin = np.zeros([len(PI)],dtype=float)
        for i in range(len(PI)):
            if Degrad == 'IZ':
                gam, Gratio[i], damping[i] = dc.IZ_1993(gam = gam,
                                                        PI=PI[i],
                                                        sigm=sigmeff[i])
            elif Degrad == 'DS':
                gam,Gratio[i],damping[i],Dmin[i] = dc.DS_2001(gam = gam,
                                                      PI=PI[i],
                                                      sigm=sigmeff[i],
                                                      OCR=OCR[i],
                                                      soil=soil[i],
                                                      N=N[i],
                                                      frq=freq[i])

        return gam, Gratio, damping,Dmin


    def _calcTauGam(self, **kwargs):
        '''
        This calculate thes stress and strains for arbitrary conditions.

        Keyword Arguments:           Definitions/Default values:
          depths              array or single value of depth at which the
                                  stresses and strains should be calculated.
                                  Default is the mid-depth of each layer.

          domain             'time' [default] or 'freq'

          MaxOnly             True [default]  or False.  When True, only the
                                  absolute maximum value of each layer is
                                  returned.  This will force the domain to be
                                  'time' and the returns to be True.

          returns             True or False [default].  When True, the stresses
                                  and strains are returned.  Otherwise, they
                                  are data attributes of the object
                                  (self.taumf and self.gammf, but only when the
                                  depths are equal to the mid-depths of the
                                  layers).

          BaseCorr            True or False [default].  Applies a baseline
                                  correction to the strain in time domain by
                                  subtracting the average value.

        Example:
            tau, gam = [DynamicsObject]._calcTauGam(depths = [0.1,0.4],
                                                domain = 'time'
                                                returns = True)

        Data attributes possibly created by this method:
        self.taumf
        self.gammf
        '''
        ComIter = self.ComIter
        t_mid = self.SoilObj.t_mid
        A = self.A
        B = self.B
        ks = self.ks
        t = self.SoilObj.t
        Gstr = self.Gstr

        depths = kwargs.get('depths',t_mid)
        domain = kwargs.get('domain','time')
        MaxOnly = kwargs.get('MaxOnly',True)
        returns = kwargs.get('returns',False)
        BaseCorr = kwargs.get('BaseCorr',False)
        if MaxOnly == True:
            domain = 'time'
            returns = True

        def hardway(t):
            '''
            Gets the requested values when things are trickier.
            returns tau, gam
            '''
            returns = True
            t = np.append(t, 1000)
            gam = np.zeros([len(depths), np.shape(A)[1]], dtype='complex128')
            tau = np.zeros([len(depths), np.shape(A)[1]], dtype='complex128')
            tsum = np.cumsum(t, dtype=float)

            j = 0
            for i in range(len(depths)):
                while depths[i] >= round(tsum[j], 3):
                    j += 1
                if depths[i] < t[0]:
                    d = depths[i]
                else:
                    d = depths[i] - tsum[j-1]
                gam[i] = 1j * ks[j] * (A[j, :] * np.exp(1j *
                            ks[j] * d) - B[j, :] *
                            np.exp(-1j * ks[j] * d))
                tau[i] = gam[i] * Gstr[j]
            return tau,gam


        # Check if I have to do it the hard way.
        if len(np.ravel(depths)) != len(np.ravel(t_mid)):
            tau,gam = hardway(t)
        else:
            if np.equal(depths, t_mid).all():
                if ComIter == True:  # If the iterations are finished, I can call existing taumf and gammf
                    gam = self.gammf
                    tau = self.taumf
                else:
                    gam = (1j * ks[:-1] *
                          (A[:-1,:] * np.exp(1j * ks[:-1] *
                           t[:,np.newaxis] /2) -
                           B[:-1,:] * np.exp(-1j * ks[:-1] *
                           t[:,np.newaxis] /2)))
                    tau = gam * Gstr[:-1,np.newaxis]
                    self.gammf = gam
                    self.taumf = tau
            else:
                tau,gam = hardway(t)

        # Do I need to convert my motions to the time domain?
        if domain == 'time':
            NPTS = self.MotObj.NPTS
            tau = np.fft.irfft(tau, axis=1)
            gam = np.fft.irfft(gam, axis=1)
            tau = tau[:,:NPTS]
            gam = gam[:,:NPTS]
            if BaseCorr == True:  # This correction is pretty elementary; I don't like it.
                basecorrect = np.mean(gam, axis=1)
                gam = gam - basecorrect[:,np.newaxis]  # Baseline correction

        if MaxOnly == True:
            tau = abs(tau).max(1)
            gam = abs(gam).max(1)
        if returns == True:
            return tau, gam


    def _AccVelDisp(self,**kwargs):
        '''
        Calculates acceleration, velocity, or displacement.

        Keyword argument:     Options:  (Default listed first)
            which          =   'Acc'  ,  'Vel'   or 'Disp'

            domain         =   'time'   or   'freq'

            MotType        =   'within', 'outcrop'  or  'incoming'

            MaxOnly        =   False  or True

            depths         =   mid-depths of each layer or list of desired

        It may be easier to use the self.calcAcc, self.calcDisp, or
            self.calcVel methods.  This is meant for internal use only.
        '''
        Which = kwargs.get('which','Acc')
        domain = kwargs.get('domain','time')
        MotType = kwargs.get('MotType','within')
        MaxOnly = kwargs.get('MaxOnly',False)
        depths = kwargs.get('depths',self.SoilObj.t_mid)
        smooth = False
        if domain == 'freq':
            smooth = kwargs.get('smooth',False)
            window = kwargs.get('window',20)
        A = self.A
        B = self.B
        h = self.SoilObj.t
        h = np.append(h,1000)
        w = self.MotObj.w
        NPTS = self.MotObj.NPTS
        g = self.g
        ks = self.ks
        disp = np.zeros([len(depths),np.shape(A)[1]], dtype='complex128')
        hsum = np.cumsum(h, dtype=float)

        # Figure out which layer we are in, and calc the displacement
        j = 0
        for i in range(len(depths)):
            while depths[i] > round(hsum[j],3):
                j += 1
            if depths[i] < h[0]:
                d = depths[i]
            else:
                d = depths[i] - hsum[j-1]
            if MotType == 'incoming':
                disp[i] = A[j,:] * np.exp(1j * ks[j] * d)
            elif MotType == 'within':
                disp[i] = (A[j,:] * np.exp(1j * ks[j] * d) +
                        B[j,:] * np.exp(-1j * ks[j] * d))
            else:
                disp[i] = 2 * A[j,:] * np.exp(1j * ks[j] * d)

        # Convert to Acc or Vel, if needed.
        if Which == 'Acc':
            out = disp * (1j * w[np.newaxis,:]) ** 2 / g
        elif Which == 'Vel':
            out = disp * (1j * w[np.newaxis,:])
        else:
            out = disp


        if domain == 'time':
            out = np.fft.irfft(out, axis=1)
            out = np.delete(out, np.s_[NPTS:], axis=1)

        # Not Fully implemented.
        if smooth == True:
            del i, j

            ind = np.arange(out.shape[1],dtype=int)
            indlow = ind - window / 2
            indhigh = ind + window / 2
            np.putmask(indlow,indlow<0,0)
            np.putmask(indhigh, indhigh>ind[-1],ind[-1])
            smoothout = np.zeros(out.shape,dtype=complex)
            for j in ind:
                smoothout[:,j] = np.mean(out[:,indlow[j]:indhigh[j]]) * (window)
            out = smoothout

        if MaxOnly == True:
            out = abs(out).max(1)

        return out


    def calcStress(self, *args, **kwargs):
        '''
        Use this method to obtain the stress at one or more depths.

        Keyword Argumentss:           Definitions/Default values:
          depths              array or single value of depth at which the
                                  stresses should be calculated.
                                  Default is the mid-depth of each layer.

          domain             'time' [default] or 'freq'

          MaxOnly             True or False [default] .  When True, only
                                  the absolute maximum value of each layer is
                                  returned.  This will force the domain to be
                                  'time' and the returns to be True.

         ratio               True or False [default].
                                 Calculates the cyclic stress ratio or CSR.

         ratiotype          Choose 'verteff' or 'meaneff' {not implemented yet}

         time                 True or False[default]. Returns an array of the
                                 time if in the time domain. Else it returns
                                 the frequency array if the domain is frequency.

        This method calls the self._calcTauGam method. It does not set any
        data attributes.
        '''
        domain = kwargs.get('domain','time')
        MaxOnly = kwargs.get('MaxOnly',False)
        ratio = kwargs.get('ratio',False)
        ratiotype = kwargs.get('ratiotype')
        time = kwargs.get('time',False)
        if ratiotype:
            ratio = True
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)

        out, strain = self._calcTauGam(depths=depths,
                                   domain=domain,
                                   MaxOnly=MaxOnly,
                                   returns=True)
        if ratio == True:
            if ratiotype == 'meaneff':
                pass
            else:
                effstr = self.SoilObj.calcInsituStress(depths, out='eff')

            if len(depths) > 1:
                out = out / effstr[:,np.newaxis]
            else:
                out = out / effstr

        if len(depths) == 1:
            out = np.ravel(out)
        if time == True:
            if domain == 'time':
                xarr = self.MotObj.time
            else:
                xarr = self.MotObj.f

            return xarr, out
        else:
            return out


    def calcStrain(self, *args, **kwargs):
        '''
        Use this method to obtain the strain at one or more depths.

        Keyword Argumentss:           Definitions/Default values:
          depths              array or single value of depth at which the
                                  strains should be calculated.
                                  Default is the mid-depth of each layer.

          domain             'time' [default] or 'freq'

          MaxOnly             True [default]  or False.  When True, only
                                  the absolute maximum value of each layer is
                                  returned.  This will force the domain to be
                                  'time' and the returns to be True.

         time                 True or False[default]. Returns an array of the
                                 time if in the time domain. Else it returns
                                 the frequency array if the domain is frequency.

        This method calls the self._calcTauGam method. It does not set any
        data attributes.
        '''
        domain = kwargs.get('domain','time')
        MaxOnly = kwargs.get('MaxOnly',False)
        time = kwargs.get('time',False)
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)

        stress, out = self._calcTauGam(depths=depths,
                                   domain=domain,
                                   MaxOnly=MaxOnly,
                                   returns=True)

        if len(depths) == 1:
            out = np.ravel(out)
        if time == True:
            if domain == 'time':
                xarr = self.MotObj.time
            else:
                xarr = self.MotObj.f

            return xarr, out
        else:
            return out


    def calcNeq(self,*args,**kwargs):
        '''
        Calculates the equivalent number of cycles using dissipated energy.

         Positional Arguments:
            depths at which equivalent number of cycles are calculated.
            ** Before refactoring, I allowed the args to be the strain ratio.
            Hopefully I have caught and changed all the old ones.

        Keyword Arguments          Definition/Default
         depths    =      array of depths, default is mid-depth of each layer
                            If the depths are specified as arguments, these
                            depths are not used.

         set_ratio  =   ratio of uniform sinusoidal amplitude and max strain
                             Note: this sets the ratio_Neq attribute and
                             defaults to the ratio used in the iteration
                             method self.ratio.  Allows me to look at how
                             the strain ratio affects the Neq in this calc.

         returns   = False {default} or True. By default, this function does
                         not return any values. Instead, the data attributes
                         of N_eq, tau_avgNeq, Neq_ratio, and dissEn are created
                         or updated.

        Data attributes created/updated by this method if he mid-layer depths
        are used. Note that these values will change for every set_ratio value.
        self.N_eq
        self.tau_avgNeq = the amplitude of the equiv. sinusoidal stress motions
        self.Neq_ratio = the strain ratio corresponding to the self.N_eq
        self.dissEn  -> only created if the mid-layer depths are used.

        This method calls the self.calcDissEn method.
        '''
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)
        ratio = kwargs.get('set_ratio',self.ratio)
        returns = kwargs.get('returns',False)
        t_mid = self.SoilObj.t_mid
        t = self.SoilObj.t

        def hardway(t):
            '''
            Gets the stiffness and damping if the mid-layer depths aren't
            specified.
            '''
            t = np.append(t, 1000)
            G = np.zeros([len(depths),], dtype=float)
            D = np.zeros([len(depths),], dtype=float)
            tsum = np.cumsum(t, dtype=float)

            j = 0
            for i in range(len(depths)):
                while depths[i] >= round(tsum[j], 3):
                    j += 1
                G[i] = self.G[j]
                D[i] = self.D[j]
            return G,D

        # Check if I have to use the hard way to get G, D.
        if len(np.ravel(depths)) != len(np.ravel(t_mid)):
            same = False
            G,D = hardway(t)

        else:
            if np.equal(depths, t_mid).all():
                same = True
                G = self.G[:-1]
                D = self.D[:-1]
            else:
                same = False
                G,D = hardway(t)

        tau, gam = self._calcTauGam(domain='time',
                                    MaxOnly=False,
                                    depths=depths,
                                    returns = True)
        DE = self.calcDissEn(tau, gam, output='normal')


        # get max stresses
        taumax = np.max(np.abs(tau),axis=1)
        tau_avgNeq = ratio * taumax  # I'm using G* here.  Not G.  Is this okay?
        dissEn_one = 2 * np.pi * D * tau_avgNeq ** 2 / G
        N_eq = DE / dissEn_one
        if same == True:
            self.tau_avgNeq = tau_avgNeq
            self.N_eq = N_eq
            self.Neq_ratio = ratio
            self.dissEn = DE

        if returns == True:
            return N_eq, tau_avgNeq, ratio, DE


    def calcDissEn(self, stresses, strains, **kwargs):
        '''
        Calculates the dissipated energy for a given set of stress and strain
        time histories.


        Input:
            stresses   m x n array of stress time histories with each layer
                            in a different row.  Make sure to use
                            stress = ifft(G* * strain(omega))
                            where G* is the complex shear modulus and
                            strain(omega) is the shear strain as a decimal in
                            the frequency domain. m is the number of layers and
                            n is the number of pts in the time history
                            (self.MotObj.NPTS)

            strains    m x n array of strain time histories (decimal),
                            in the same format as the stresses.
        Keyword Arguments:

            FinalValOnly = True | False


        Output:
            DissE    The dissipated energy in units of stress per unit volume
                        if output='normal'.  Else, an array of strain ratios.
                        To obtain the normalized dissipated energy like that
                        given by SHAKEVT, divide the dissipated energy
                        by the vertical effective stress.
        '''
        if kwargs.get('FinalValOnly', True):
            dissEn = np.zeros(np.shape(stresses)[0], dtype=float)
            dissEn =  np.sum((stresses[:,1:] + stresses[:,:-1]) *
                            (strains[:,1:] - strains[:,:-1]), axis = 1) * 0.5
            return dissEn
        else:
            dissEn = np.zeros(stresses.shape, dtype=float)
            dissEn =  np.cumsum((stresses[:,1:] + stresses[:,:-1]) *
                                (strains[:,1:] - strains[:,:-1]) * 0.5,
                                axis = 1)
            return dissEn



    def calcStrainRatio_DE(self, stresses, strains, **kwargs):
        '''
        Calculates a strain ratio using the weighted average
                        strain from the dissipated energy.

        Input:
            stresses   m x n array of stress time histories with each layer
                            in a different row.  Make sure to use
                            stress = ifft(G* * strain(omega))
                            where G* is the complex shear modulus and
                            strain(omega) is the shear strain as a decimal in
                            the frequency domain. m is the number of layers and
                            n is the number of pts in the time history
                            (self.MotObj.NPTS)

            strains    m x n array of strain time histories (decimal),
                            in the same format as the stresses.

        Output: An array of strain ratios of length m.

        '''
        plot = kwargs.get('plot', False)
        dissEn = np.zeros(np.shape(stresses)[0], dtype=float)
        dissEn =  ((stresses[:,1:] + stresses[:,:-1]) *
    					(strains[:,1:] - strains[:,:-1])) * 0.5
        x = stresses/strains
        dx = np.diff(x)
        avgstrain = np.zeros(strains.shape[0],)

        for i in range(np.shape(x)[0]):
            ind = np.where((dx[i,:-2] < 0) & (dx[i,1:-1] > 0) & (dx[i,2:]<0))
            if ind[0][0] == 0:
                ind = ind[0] + 1
            else:
                ind = np.insert(ind[0],0,0)
            dDE = np.zeros((ind.shape[0]),)
            dstrain = np.zeros((ind.shape[0]),)
#            dstress = np.zeros((ind.shape[0]),)

            for k in range(len(ind)-2):
                if k % 2 == 1:
                    continue
                else:
                    dstrain[k] = (np.max(strains[i,ind[k]:ind[k+2]]) -
                                  np.min(strains[i,ind[k]:ind[k+2]]))/2
#                    dstress[k] = (np.max(stresses[i,ind[k]:ind[k+2]]) -
#                                  np.min(stresses[i,ind[k]:ind[k+2]]))/2
                    dDE[k] = dissEn[i,ind[k+2]] - dissEn[i,ind[k]]

                    if plot == True:
                        plt.plot(strains[i,ind[k]:ind[k+2]],
                                 stresses[i,ind[k]:ind[k+2]],'-b')
                        plt.xlabel('Strain')
                        plt.ylabel('Stress (kPa)')
                        plt.show()
            avgstrain[i] = np.average(dstrain,weights=dDE)
#            G = dstress / dstrain
#            print(self.G[i])
#            pdb.set_trace()
            del dDE, ind, dstrain#, dstress
        strainratio = avgstrain / np.max(strains,axis=1)
        np.putmask(strainratio,strainratio<=0,0.1)
        np.putmask(strainratio,strainratio>1.,1.)
        return strainratio


    def calcDisp(self,*args,**kwargs):
        '''
        Calculates displacements using the _AccVelDisp method.

        Positional Arguments, if used, are the depths.

        Keyword argument:     Options:  (Default listed first)

            domain         =   'time'   or   'freq'

            MotType        =   'within', 'outcrop'  or  'incoming'

            MaxOnly        =   False  or True

            depths         =   mid-depths of each layer or list of desired
                                These depths can also be listed as regular,
                                non-keyword arguments; see the example.

         Returns an array of velocities in units of m for time
                 histories, m-s for frequency domain. Array is 2D if more
                 than one depth is specified.

         Example:
             disp = DynamicObject.calcDisp(2.0,4.5,
                                       domain='time',
                                       MaxOnly=False,
                                       MotType = 'within')
        '''
        domain = kwargs.get('domain','time')
        MotType = kwargs.get('MotType','within')
        MaxOnly = kwargs.get('MaxOnly',False)

        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)

        out = self._AccVelDisp(which='Disp',depths=depths,domain=domain,
                                MotType=MotType, MaxOnly=MaxOnly)
        if len(depths) == 1:
            out = np.ravel(out)
        return out


    def calcVel(self,*args,**kwargs):
        '''
        Calculates velocity using the _AccVelDisp method.

        Keyword argument:     Options:  (Default listed first)

            domain         =   'time'   or   'freq'

            MotType        =   'within', 'outcrop'  or  'incoming'

            MaxOnly        =   False  or True

            depths         =   mid-depths of each layer or list of desired
                                 These depths can also be listed as regular,
                                 non-keyword arguments; see the example.

             Returns an array of velocities in units of m/s for time
                 histories, m/s-s for frequency domain. Array is 2D if more
                 than one depth is specified.


         Example:
             vel = DynamicObject.calcVel(2.0,4.5,
                                       domain='time',
                                       MaxOnly=False,
                                       MotType = 'within')


        '''
        domain = kwargs.get('domain','time')
        MotType = kwargs.get('MotType','within')
        MaxOnly = kwargs.get('MaxOnly',False)
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)

        out = self._AccVelDisp(which='Vel',depths=depths,domain=domain,
                                MotType=MotType, MaxOnly=MaxOnly)
        if len(depths) == 1:
            out = np.ravel(out)
        return out


    def calcAcc(self,*args,**kwargs):
        '''
        Calculates acceleration using the _AccVelDisp method.

        Keyword argument:     Options:  (Default listed first)

            domain         =   'time'   or   'freq'

            MotType        =   'within', 'outcrop'  or  'incoming'

            MaxOnly        =   False  or True

            depths         =   mid-depths of each layer or list of desired
                                 These depths can also be listed as regular,
                                 non-keyword arguments; see the example.

         Returns an array of accelerations.

         Example:
             acc = DynamicObject.calcAcc(2.0,4.5,
                                       domain='time',
                                       MaxOnly=False,
                                       MotType = 'within')

            Returns an array of velocities in units of m/s/s for time
                 histories, m/s/s-s for frequency domain. Array is 2D if more
                 than one depth is specified.
        '''
        domain = kwargs.get('domain','time')
        MotType = kwargs.get('MotType','within')
        MaxOnly = kwargs.get('MaxOnly',False)
        smooth = kwargs.get('smooth',False)
        window = kwargs.get('window',20.)
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)

        out = self._AccVelDisp(which='Acc',depths=depths,domain=domain,
                                MotType=MotType, MaxOnly=MaxOnly, smooth=smooth,
                                window=window)
        if len(depths) == 1:
            out = np.ravel(out)
        return out


    def TF(self,*args, **kwargs):
        '''
        Calculates the transfer function between two layers in the frequency
            domain.

        Positional Arguments:
            List  of two layers. A 1 refers to the top layer.
                                    0 refers to bedrock.
        Keyword Arguments:
            MotType =  'outcrop'[default] or 'within',
                        'rocktosurf','rocktobase'
        Output
            tf  = the transfer function.
        '''
        args = list(args)
        layer1 = args[0]
        layer2 = args[1]
        MotType = kwargs.get('MotType','outcrop')
        A = self.A
        B = self.B
        layer1 -= 1
        layer2 -= 1

        if MotType == 'rocktosurf':
            tf = np.abs((A[0,:] + B[0,:])/(2 * A[-1,:]))
        elif MotType == 'rocktobase':
            tf = np.abs((A[-1,:] + B[-1,:])/(2 * A[-1,:]))
        elif MotType == 'outcrop':
            tf = np.abs(A[layer1,:]/A[layer2,:])
        else:
            tf = np.abs((A[layer1,:]+B[layer1,:])/
                        (A[layer2,:]+B[layer2,:]))
        return tf



    def calcTF(self,*args, **kwargs):
        '''
        Calculates the transfer function between two layers in the frequency
            domain.

        Positional Arguments:
            List  of two layers. A 1 refers to the top layer.
                                    0 refers to bedrock.
        Keyword Arguments:
            MotType =  Tuple of size two with values of
                            'outcrop' or 'within', [default]

        Output
            tf  = the transfer function.

        Example:
            tf = self.calcTF(1,0,MotType=('within', 'outcrop'))
        '''
        args = list(args)
        layer1 = args[0]
        layer2 = args[1]
        MotType = kwargs.get('MotType','outcrop')
        A = self.A
        B = self.B
        layer1 -= 1
        layer2 -= 1

        if MotType[0] == 'outcrop':
            num = 2 * A[layer1,:]
        else:
            num = A[layer1,:] + B[layer1,:]

        if MotType[1] == 'outcrop':
            den = 2 * A[layer2,:]
        else:
            den = A[layer2,:] + B[layer2,:]

        return np.abs(num/den)



    def _PiecewiseExact(self,acc_g,w,dt,zeta,**kwargs):
        '''
        Determines the spectral response of an SDOF structure.

        Uses Duhamel's integral.  *** This is currently really slow!!!***

        Positional Arguments:
            acc_g = ground acceleration in units of g

            w = natural frequency of the SDOF structure.
                To return an acceleration response vector,
                input a vector of frequencies here.

            dt = sampling interval or time interval of the input acceleration

            zeta = the damping ratio of the structure = c/c_crit (decimal)

        Keyword Arguments:

            returnextra     If True, returns displacement, velocity, and
                                total acceleration. Only useful for single
                                value of w. Else, returns Spectral
                                Acceleration.  Default is False.

        Code from rsc.c at
        <http://pubs.usgs.gov/of/2006/1369/ResponseSpectra.zip>
        This would be a good function to move to C.
        '''
        returnextra = kwargs.get('returnextra', False)
        output_name = kwargs.get('output', 'PSa')
        a = np.float64(np.ravel(acc_g) * 9.80665) # m/s/s
        w = np.float64(np.ravel(w))
        dt = np.float64(dt)
        z = np.float64(zeta)


        root = np.sqrt(1. - z ** 2)
        wd = w * root
        sine = np.sin(wd * dt)
        cosn = np.cos(wd * dt)
        expo = np.exp(-z * w * dt)

        zratio = z / root
        blue = 2. * z / (w ** 3 * dt)
        yell = (2 * z ** 2 - 1) / (w ** 2 * dt)

        a11 = expo * (zratio * sine + cosn)
        a12 = expo * sine / wd
        a21 = -w / root * expo * sine
        a22 = expo * (cosn - zratio * sine)

        b11 = (expo * ((yell + z / w) * (sine / wd)
				+ (blue + 1. / w ** 2) * cosn) - blue)
        b12 = (-expo * (yell * sine / wd + blue * cosn)
				- (1. / w ** 2) +	blue)
        b21 = (expo * ((yell + z / w) * (cosn - zratio * sine)
				-(blue + 1. / w ** 2) * (wd * sine + z * w * cosn)) +
				1. / (w ** 2 * dt))
        b22 = (-expo * (yell * (cosn - zratio * sine) -
				blue * (wd * sine + z * w * cosn)) -
				1. / (w ** 2 * dt))
        if np.isscalar(w):
            u = np.zeros([1,len(a)], dtype=np.float64)
            v = np.zeros([1,len(a)], dtype=np.float64)
            for i in range(len(a) - 1):
                u[:,i+1] = (a11 * u[:,i] + a12 * v[:,i] + b11 * a[i] +
                            b12 * a[i+1])
                v[:,i+1] = (a21 * u[:,i] + a22 * v[:,i] + b21 * a[i] +
                            b22 * a[i+1])
            if output_name is 'Sa':                
                output = - (2 * z * w * v + w ** 2 * u) / 9.80665
            elif output_name is 'PSa':
                output = - ( w ** 2 * u) / 9.80665
            # elif output_name is 'Sv':
            #     output = - ([0,u[1:]-u[:-1]]/dt)
            elif output_name is 'PSv':
                output = - ( w * u)
            elif output_name is 'Sd':
                output = - (u) 

                
        else:
            u = np.zeros([len(w),len(a)], dtype=np.float64)
            v = np.zeros([len(w),len(a)], dtype=np.float64)
            for i in range(len(a) - 1):
                u[:,i+1] = (a11 * u[:,i] + a12 * v[:,i] + b11 * a[i] +
                            b12 * a[i+1])
                v[:,i+1] = (a21 * u[:,i] + a22 * v[:,i] + b21 * a[i] +
                            b22 * a[i+1])
            if output_name is 'Sa':                
                output = - (2 * z * w[:,np.newaxis] * v + w[:,np.newaxis] ** 2 * u) / 9.80665
            elif output_name is 'PSa':
                output = - ( w[:,np.newaxis] ** 2 * u) / 9.80665
            # elif output_name is 'Sv':
            #     output = - ([0,u[1:]-u[:-1]]/dt)
            elif output_name is 'PSv':
                output = - ( w[:,np.newaxis] * u)
            elif output_name is 'Sd':
                output = - (u) 

        if returnextra:
            return np.ravel(u), np.ravel(v), np.ravel(TA)
        else:
            output = np.max(np.abs(output),axis=1)
            return output
        
    def fouriersolution(self,FA_g,f,w,zeta,**kwargs):
        
        FA = FA_g * 9.80665 # m/s/s
        w = np.float64(np.ravel(w))
        z = np.float64(zeta)
        ii=0
        Sa = np.ones(len(w))
        for w1 in w:        
            fn = w1/(np.pi*2)
            hf = -fn**2/(f**2-fn**2-2*1j*z*f*fn)
            acc_fr_s = FA*hf
            acc_t_s = np.fft.irfft(acc_fr_s)
            Sa[ii] = np.max(np.absolute(acc_t_s))
            ii=ii+1
        return Sa    


    def StressRedCoef(self, *args, **kwargs):
        '''
        Deprecated name.  Use calc_r_d instead.
        '''
        print('Function name deprecated.' +
                '    Use "calc_r_d" instead of "StressRedCoef"')
        if kwargs.get('returns', False):
            rd = self.calc_r_d(*args, **kwargs)
            return rd
        else:
            self.calc_r_d(*args, **kwargs)


    def calc_r_d(self, *args, **kwargs):
        '''
        Calculates the stress reduction coefficient for the profile and
        ground motion.

        Positional Arguments:
        Depths at which the r_d coefficient should be calculated.

        Keyword Arguments:          Definition/Default
        depths                      Default is the mid-depth of each layer.
                                        Disregarded if positional arguments
                                        are used.

        returns                     False [default] or True.

        If returns == True, stress reduction coefficient (r_d) is returned.
        This will return the calculated values.

        If the mid-layer depths are used, the self.rd data attribute is created
        by this method.
        '''
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)
        returns = kwargs.get('returns',False)
        if np.equal(depths, self.SoilObj.t_mid).all():
            same = True
            sigv = self.SoilObj.sigv
        else:
            pdb.set_trace()
            same = False
            sigv = self.SoilObj.calcInsituStress(*depths,out='total')

        taumax,gammax = self._calcTauGam(depths=depths,domain='time',
                               MaxOnly=True,returns=True)
        amax = self.calcAcc(depths=[0.],MaxOnly=True,MotType='outcrop',
                        domain='time')

        rd = taumax / (amax * sigv)
        if same == True:
            self.rd = rd
        if returns == True:
            return rd


    def RespSpec(self, *args, **kwargs):
        '''
        Calculates the response spectra for a given depth(s).

        Positional Arguments:
            Depths at which RS is to be calculated. If none are given,
                it will calculate the response spectra for
                all mid-depths (very slow).

        Keyword Arguments:        Values:
        depths                an array of depths.  Disregarded if the
                                    positional arguments are used.

        damping               any damping value as decimal, default is 0.05

        w                     an array of frequencies. Computationally
                                intensive if this array is large.
        method               'piecewiseexact' is the only one implemented now.
        f                    frequency of FA

        MotType        =   'within', 'outcrop'  or  'incoming'

        Returns Sa, w

        This method calls the self._PiecewiseExact method.
        '''
        if len(args) > 0:
            depths = np.array(args)
        else:
            depths = kwargs.get('depths',self.SoilObj.t_mid)
        damping = kwargs.get('damping',0.05)
        # T="0.01  0.02    0.022   0.025   0.029   0.03    0.032   0.035   0.036   0.04    0.042   0.044   0.045   0.046   0.048   0.05    0.055   0.06    0.065   0.067   0.07    0.075   0.08    0.085   0.09    0.095   0.1 0.11    0.12    0.13    0.133   0.14    0.15    0.16    0.17    0.18    0.19    0.2 0.22    0.24    0.25    0.26    0.28    0.29    0.3 0.32    0.34    0.35    0.36    0.38    0.4 0.42    0.44    0.45    0.46    0.48    0.5 0.55    0.6 0.65    0.667   0.7 0.75    0.8 0.85    0.9 0.95    1   1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2   2.2 2.4 2.5 2.6 2.8 3   3.2 3.4 3.5 3.6 3.8 4   4.2 4.4 4.6 4.8 5   5.5 6   6.5 7   7.5 8   8.5 9   9.5 10  11  12  13  14  15  20"
        # T=[float(i) for i in T.split()]
        # T = kwargs.get('T', np.arange(T))
        # f = 1/T
        # w = 2 * np.pi * f
        w = kwargs.get('w',np.arange(0.001,600,1.2))
        f = kwargs.get('f',np.arange(0.001,600,1.2))#<---- 
        method = kwargs.get('method', 'piecewiseexact')
        MotType = kwargs.get('MotType','within')

        if method == 'piecewiseexact':
            acc = self.calcAcc(depths=depths,
                           domain='time',
                           MaxOnly=False,
                           MotType = MotType)
        

            SA = np.zeros([len(depths),len(w)],dtype=float)
            for i in range(len(depths)):
                if len(depths) == 1:
                    SA = self._PiecewiseExact(acc, w, self.MotObj.dt,damping)
                else:
                    SA[i,:] = self._PiecewiseExact(acc[i,:], w,
                                                self.MotObj.dt,damping)
            if SA.shape[0] == 1:
                SA = np.ravel(SA)
            return SA,w
            # return SA,f
        elif method == 'fouriersolution':
            FA = self.calcAcc(depths=depths,
               domain='freq',
               MaxOnly=False,
               MotType = MotType)
            SA = np.zeros([len(depths),len(w)],dtype=float)
            for i in range(len(depths)):
                if len(depths) == 1:
                    SA = self.fouriersolution(FA, w,damping)
                else:
                    SA[i,:] = self.fouriersolution(FA[i,:], w, damping)
            if SA.shape[0] == 1:
                SA = np.ravel(SA)
            return SA/9.8,w
            # return SA/9.8,f
        else:
            print('Did not understand the method. Please try again.')
            return
        


    def ASCE_RS_ratios(self,**kwargs):
        '''
        Calculates the ratio between the 5% damped surface
        response spectra of surface ground motions to input base motions.

        See ASCE 7-05, sec. 21.1.3.

        Keyword Arguments:
            T       array of periods

        Returns the surface to base response ratio and the corresponding
            array of periods.

        This method calls the self.RespSpec method and assumes 5% damping,
            'within' surface motion, and 'outcrop' base motion.

        Example:

            RS_ratio, T = self.ASCE_RS_ratios(T=np.logspace(-3,1,num=1000))
        '''
        T = kwargs.get('T',np.logspace(-3,1,num=1000))

        # Obtain the surface response spectrum
        surf,w = self.RespSpec(0.,damping=0.05,w=2.*np.pi/T,MotType='within')

        # Obtain the input base motion response spectrum
        basedepth = np.sum(self.SoilObj.t) + 0.001
        base,w = self.RespSpec(basedepth,damping=0.05,w=w,MotType='outcrop')

        ratio = surf / base

        return ratio,T


    def plot(self, *args, **kwargs):
        '''
        Outputs standard plots that illustrate the equivalent linear response.

        Positional Arguments:
            List the plots you want to see.
            Possible plots:
                mod         modulus profile
                damp        damping profile
                maxacc      max acceleration profile
                csr         max stress ratio profile
                maxstrain   max strain profile
                dissen      dissipated energy profile
                                (calc'd only at mid points)
                strainratio strain ratio profile used to obtain average strain
                                from max
                cycles      number of equivalent cycles profile
                rd          stress reduction coefficient
                allprofile  or all   outputs all



        For keywords arguments:
            dz =    For profiles: depth increments (in meters)

            acc =   an array containing the following:
                        [depth,'time' or 'freq','outcrop' or 'within']

            vel =   an array containing the following:
                        [depth,'time' or 'freq','outcrop' or 'within']

            disp =  an array containing the following:
                        [depth,'time' or 'freq','outcrop' or 'within']

            tf   =  an array containing the following:
                        [toplayer, bottomlayer, 'outcrop' or 'within']
                        Use the layer numbers.  If only one value is given,
                        the bottom layer is assumed to be bedrock.

            respspec = an array containing the following:
                        [depth,damping,'outcrop' or 'within']

            plottype = 'loglog' [default], 'semilogx',
                            'linear' (only for freq domain)

            period = False [default] or True. Specifies if the horizontal axis
                     should show period or frequency.

            save = False (default) or True. For now, if save = True, then the
                    figures will also be closed. (Haven't figured it out yet.)

            extension = '.png' (default) or any other that matplotlib handles

            syms1   =  matplotlib line symbols.  Default is '-b'

            syms2   =  matplotlib line symbols.  Default is 'ob'

            linelabel =   label for plot.  Default is ''


        Example:

            DynamicsObj.plot('mod','maxacc', save=True, dz=1)
        '''
        syms1 = kwargs.get('syms1','-b')
        syms2 = kwargs.get('syms2','ob')
        linelabel = kwargs.get('linelabel','')
        t = self.SoilObj.t
        t_mid = self.SoilObj.t_mid
        basedepth = np.sum(t,dtype=float)

        if 'all'  in args:
            args = ('mod','damp','maxacc','csr','maxstrain',
                        'dissen','strainratio','cycles', 'rd')
        if 'allprofile' in args:
            args = ('mod','damp','maxacc','csr','maxstrain',
                        'dissen','strainratio','cycles', 'rd')

        dz = kwargs.get('dz')
        acc = kwargs.get('acc')
        vel = kwargs.get('vel')
        disp = kwargs.get('disp')
        tf = kwargs.get('tf')
        respspec = kwargs.get('respspec')
        period = kwargs.get('period',False)
        plottype = kwargs.get('plottype','loglog')
        save = kwargs.get('save',False)
        ext = kwargs.get('extension','.png')
        return_fig = kwargs.get('fig')
        if return_fig != None:
            plt.figure(return_fig.number)

        if dz:
            depths = np.arange(0.,basedepth+dz,dz,dtype=float)
        else:
            depths = t_mid


        def profilecreator(array, **kwargs):
            '''
            Makes points at the top and bottom of each layer to make plotting
            better.
            '''
            depth = kwargs.get('depth',False)
            if depth == True:
                output = [0]
                for i in range(len(array)):
                    if output[i] == 0:
                        output.append(array[i])
                    else:
                        output.append(array[i-1])
                        output.append(array[i])
                output.append(array[-1])
                output.append(array[-1]*1.05)
            else:
                output = []
                for i in range(len(array)):
                    output.append(array[i])
                    output.append(array[i])
            return np.array(output)


        def profinterp(depth2,values2,newdepths):
            newvalues = []
            for i in range(len(newdepths)):
                for j in range(1,len(depth2),2):
                    if (newdepths[i] <= depth2[j])&(newdepths[i] >= depth2[j-1]):
                        newvalues.append(values2[j])
                        break
            return newvalues


        depth2 = profilecreator(np.cumsum(t),depth=True)

        if 'mod' in args:
            mod = profilecreator(self.G)
            plt.figure('Shear Modulus Profile')
            plt.plot(mod,depth2,syms1)
            plt.plot(self.G[:-1],t_mid, syms2)
            plt.xlabel(r'Degraded Shear Modulus, $G$  $(kPa)$')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, np.max(self.G[:-1])*1.1, ymax, ymin])
            fig = plt.gcf()

        if 'damp' in args:
            damp = profilecreator(self.D)
            plt.figure('Damping Profile')
            plt.plot(damp*100,depth2, syms1)
            plt.plot(self.D[:-1]*100,t_mid,'or')
            plt.xlabel(r'Degraded Damping, $D$  ($\%$)')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()

        if 'dissen' in args:
            if self.dissEn != []:
                dissen = self.dissEn
            else:
                self.calcNeq()
                dissen = self.dissEn
            dissen2 = profilecreator(dissen)
            plt.figure('DissEn')
            plt.plot(dissen2,depth2[:-2], syms1)
            plt.plot(self.dissEn,t_mid,'or')
            plt.xlabel(r'Dissipated Energy  $(kJ/m^3)$')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()

        if 'cycles' in args:
            try:
                cyc = self.N_eq
            except AttributeError:
                self.calcNeq()
                cyc = self.N_eq
            cyc2 = profilecreator(cyc)
            plt.figure('NumCyc')
            plt.plot(cyc2,depth2[:-2], syms1)
            plt.plot(self.N_eq,t_mid,'or')
            plt.xlabel(
                'Number of Equivalent Cycles Based on Dissipated Energy')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()

        if 'maxacc' in args:
            accel = self.calcAcc(domain='time',
                             MotType='within',
                             MaxOnly=True,
                             depths=depths)
            plt.figure('AccProfile')
            plt.plot(accel,depths, syms1)
            plt.xlabel(r'Max Acceleration ($g$)')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()

        if ('maxstrain' in args) or ('csr' in args):
            tau,gam = self._calcTauGam(domain='time',
                             returns=True,
                             MaxOnly=True,
                             depths=depths)

        if 'maxstrain' in args:
            plt.figure('MaxStrain')
            plt.plot(gam*100,depths, syms1)
            plt.xlabel(r'Max Shear Strain, $\gamma$, ($\%$)')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()

        if 'csr' in args:
            Gnew = np.array(profinterp(depth2,mod,depths))
            sigveff = self.SoilObj.calcInsituStress(depths,out='eff')
            CSR = gam * Gnew / sigveff
            plt.figure('CSR')
            plt.plot(CSR,depths, syms1)
            plt.xlabel(r'Max Cyclic Stress Ratio, $CSR$')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()

        if 'strainratio' in args:
            ratio1 = np.ones(t.shape)*np.array(self.ratio)
            ratio = profilecreator(ratio1)
            plt.figure('StrainRatio')
            plt.plot(ratio,depth2[:-2], syms1)
#            plt.plot(ratio1,t_mid,'or')
            plt.xlabel('Strain Ratio')
            plt.ylabel(r'Depth ($m$)')
            plt.title(self.SoilObj.filename+'-'+self.MotObj.filename)
            [xmin, xmax, ymin, ymax] = plt.axis()
            plt.axis([xmin, xmax, ymax, ymin])
            fig = plt.gcf()


        # Specific Depth  plots
        def plothelpers(x,y,domain,plottype,period,figure,dt=True,label=None):
            '''
             x should be frequency or time, not angular frequency
             If dt is True, the response is multiplied by dt.
            '''
            plt.figure(figure.number)
            if domain == 'freq':
                if dt:
                    y = np.abs(y) * self.MotObj.dt

                if period:
                    x = 1./x
                    plt.xlabel('Period (s)')
                    if plottype == 'loglog':
                        plt.loglog(x,y,label=label)

                    elif plottype == 'semilogx':
                        plt.semilogx(x,y,label=label)
                    else:
                        plt.plot(x,y,label=label)
                else:
                    plt.xlabel(r'Frequency ($Hz$)')
                    if plottype == 'loglog':
                        plt.loglog(x,y,label=label)
                    elif plottype == 'semilogx':
                        plt.semilogx(x,y,label=label)
                    else:
                        plt.plot(x,y,label=label)
            else:
                plt.plot(x,y,label=label)
                plt.xlabel('Time (s)')
            curr_fig = plt.gcf()
            return curr_fig

        if acc:
            if acc[2] == 'outcrop':
                mot = 'outcrop'
            else:
                mot = 'within'
            Accel = np.ravel(self.calcAcc(acc[0],domain=acc[1],MotType=mot))
            if acc[1] == 'time':
                x = self.MotObj.time
            else:
                x = self.MotObj.f

            if return_fig:
                plt.figure(return_fig.number)
            else:
                plt.figure(r'Acceleration at z = %s'%str(acc[0]))
            if acc[1] == 'time':
                plt.ylabel(r'Acceleration ($g$)')
            else:
                plt.ylabel(r'Acceleration ($\frac{m}{s^2}$ $-s$)'
                            %str(acc[0]) )
            plt.title(self.SoilObj.filename+' - '+self.MotObj.filename)
            fig = plt.gcf()
            label = '%s m'%str(acc[0])
            fig = plothelpers(x,Accel,acc[1],plottype,period,fig,label=label)

        if vel:
            if vel[2] == 'outcrop':
                mot = 'outcrop'
            else:
                mot = 'within'
            Veloc = np.ravel(self.calcVel(vel[0],domain=vel[1],MotType=mot))
            if vel[1] == 'time':
                x = self.MotObj.time
            else:
                x = self.MotObj.f

            if return_fig:
                plt.figure(return_fig.number)
            else:
                plt.figure(r'Velocity at z = %s'%str(vel[0]))
            if vel[1] == 'time':
                plt.ylabel(r'Velocity ($\frac{m}{s}$)')
            else:
                plt.ylabel(r'Velocity ($\frac{m}{s}$ $-s$)')
            plt.title(self.SoilObj.filename+' - '+self.MotObj.filename)
            fig = plt.gcf()
            label = '%s m'%str(vel[0])
            fig = plothelpers(x,Veloc,vel[1],plottype,period,fig,label=label)

        if disp:
            if disp[2] == 'outcrop':
                mot = 'outcrop'
            else:
                mot = 'within'
            Disp = np.ravel(self.calcDisp(disp[0],domain=disp[1],MotType=mot))
            if disp[1] == 'time':
                x = self.MotObj.time
            else:
                x = self.MotObj.f

            if return_fig:
                plt.figure(return_fig.number)
            else:
                plt.figure(r'Displacement at z = %s'%str(disp[0]))
            if disp[1] == 'time':
                plt.ylabel(r'Displacement ($m$)')
            else:
                plt.ylabel(r'Displacement ($m$ $-s$)'
                            %str(disp[0]) )
            plt.title(self.SoilObj.filename+' - '+self.MotObj.filename)
            fig = plt.gcf()
            label = '%s m'%str(disp[0])
            fig = plothelpers(x,Disp,disp[1],plottype,period,fig,label=label)


        if tf:
            x = self.MotObj.f
            if len(tf) == 1:
                try:
                    layer1 = int(tf[0])
                    layer2 = 'Bedrock'
                except:
                    print("I didn't understand the transfer function input")
#                    break
                y = self.TF(layer1,0)
            elif len(tf) == 2:
                try:
                    layer1 = int(tf[0])
                    layer2 = int(tf[1])
                    y = self.TF(layer1,layer2)
                except:
                    print ("Assuming the second value of the TF " +
                            "specifies motion type")
                    mot = tf[1]
                    layer2 = 'Bedrock'
                    y = self.TF(layer1,0,MotType=mot)
            else:
                layer1 = int(tf[0])
                layer2 = int(tf[1])
                mot = tf[2]
                y = self.TF(layer1,layer2,MotType=mot)
            plt.figure('Transfer Function' )
            plt.ylabel(r'Transfer Function Between Layers %s and %s'
                        %(str(layer1),str(layer2)))
            plt.title(self.SoilObj.filename+' - '+self.MotObj.filename)
            fig = plt.gcf()
            label = 'Layer #%s'%str(layer1)
            fig = plothelpers(x,y,'freq',plottype,period,fig,dt=False,
                              label=label)


        if respspec:
            plt.figure('Response Spectrum' )
            plt.ylabel(r'Response Spectrum at a Depth of %s m'
                        %str(respspec[0]))
            plt.title(self.SoilObj.filename+' - '+self.MotObj.filename)
            y,w = self.RespSpec(respspec[0],
                                damping=respspec[1],
                                MotType=respspec[2],
                                w=self.MotObj.w[1:])
            x = self.MotObj.f[1:]
            y = np.ravel(y)
            fig = plt.gcf()
            label = 'Depth: %s m, Damping: %s'%(str(respspec[0]),str(respspec[1]))
            fig = plothelpers(x,y,'freq',plottype,period,fig,dt=False,
                              label=label)

        if save == True:
            import matplotlib
            fig=[manager.canvas.figure for manager in
                        matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
            for i, figure in enumerate(fig):
                figure.savefig(self.SoilObj.profile_name[:5] + '_' +
                                self.MotObj.filename[:5]+'-'+
                                figure._label+ext)
        return fig




class DSvariables(object):
    '''
    A simple object to hold data for Darendeli and Stokoe's degradation
    curve options.
    '''
    def __init__(self,kwargs,SoilObj):
        # Degradation Curve Overrides and Special Values:
        self.sigmeff = kwargs.get('DCO_sigmeff',SoilObj.sigmeff)
        self.PI = kwargs.get('DCO_PI',SoilObj.PI)
        PI = self.PI
        if len(SoilObj.OCR) == 0:
            self.OCR = kwargs.get('DS_OCR',
                                  1. * np.ones(len(PI),dtype=float)) #Non-over consolidated as default
        else:
            self.OCR = SoilObj.OCR
        self.soiltype = kwargs.get('DS_soil', SoilObj.DS_soil)
        if self.soiltype == None:
            self.soiltype = 0. * np.ones(len(PI),dtype=float)  #General default
        self.N = kwargs.get('DS_N', SoilObj.DS_N)
        if self.N == None:
            self.N = 10. * np.ones(len(PI),dtype=float) # 10 cycles default
        self.freq = kwargs.get('DS_freq', SoilObj.DS_freq)
        if self.freq == None:
            self.freq = 1. * np.ones(len(PI),dtype=float)  # 1 Hz default

def nextpow2(i):
    '''
    Returns the the log-base-2 of the next power of 2 greater than the input
    value.

    Example:
    >>> y = nextpow2(5)
    >>> print y
    3
    >>> 2**nextpow2(5)
    8
    '''
    n = 2
    while n < i:
        n = n * 2
    return int(round(np.log2(n)))


def input_soil_profile(csvprofile):
    '''
    Creates a soil profile object using Lake's format of the csv.
    The new format is simpler and doesn't contain any options for output.

    Positional Argument:
        Name or path of a profile .csv file.
    '''
    layerno = []
    layername = []
    MRD = []
    gam_c = []
    unit_wt = []
    t = []
    ko = []
    phi = []
    Vs = []
    PI = []
    N160 = []
    OCR = []
    DS_soil = []
    DS_N = []
    DS_freq = []
    Cu=[]
    D50=[]
    Su=[]

    count = 0
    filein = open(csvprofile, "rU")
    reader = csv.reader(filein)
    for line in reader:
        count += 1
        if count == 1:
            Profile_name = str(line[1])
        elif count == 2:
            units = int(line[1])
        elif count == 3:
            gwt = float(line[1])
        elif count == 4:
            loc_mot = int(line[1])
        elif count == 5:
            outcrop_mot = int(line[1])
        elif count == 6:
            pass
        else:
            layerno.append(int(line[0]))
            layername.append(line[1])
            MRD.append(line[2])
            try:
                gam_c.append(float(line[3]))
            except (IndexError,ValueError):
                gam_c.append(0.01) # If there is no number, I assume 1.

            unit_wt.append(float(line[4]))

            Vs.append(float(line[5]))
            try:
                t.append(float(line[6]))
            except (ValueError, IndexError):
                break # Notice I kick out if the thicknesses run out.
            try:
                PI.append(float(line[7]))
            except (ValueError,IndexError):
                PI.append(0.) # If there is no number, I assume 0.
            try:
               ko.append(float(line[8]))
            except (ValueError,IndexError):
               ko.append(0.) # If there is no number, I assume 0.
            try:
                phi.append(float(line[9]))
            except (ValueError,IndexError):
                phi.append(0.) # If there is no number, I assume 0.
            try:
                N160.append(float(line[10]))
            except (ValueError,IndexError):
                N160.append(0.) # If there is no number, I assume 0.
            try:
                OCR.append(float(line[11]))
            except (IndexError,ValueError):
                OCR.append(1.) # If there is no number, I assume 1.
            try:
                DS_soil.append(int(line[12]))
            except (IndexError,ValueError):
                DS_soil.append(0) # If there is no number, I assume 0, general soil.
            try:
                DS_N.append(int(line[13]))
            except (IndexError,ValueError):
                DS_N.append(10) # If there is no number, I assume 10.
            try:
                DS_freq.append(float(line[14]))
            except (IndexError,ValueError):
                DS_freq.append(1) # If there is no number, I assume 1.
            try:
                D50.append(float(line[15]))
            except (IndexError,ValueError):
                D50.append(1) # If there is no number, I assume 1.
            try:
                Cu.append(float(line[16]))
            except (IndexError,ValueError):
                Cu.append(1) # If there is no number, I assume 1.
            try:
                Su.append(float(line[17]))
            except (IndexError,ValueError):
                Su.append(-1.0) # If there is no number, I assume -1.


    filein.close()
    # Convert to numpy arrays
    layerno = np.array(layerno)
    unit_wt = np.array(unit_wt)
    ko = np.array(ko)
    Vs = np.array(Vs)
    t = np.array(t)
    phi = np.array(phi)
    PI = np.array(PI)
    N160 = np.array(N160)
    OCR = np.array(OCR)
    DS_soil = np.array(DS_soil)
    DS_N = np.array(DS_N)
    DS_freq = np.array(DS_freq)
    D50 = np.array(D50)
    Cu = np.array(Cu)
    gam_c = np.array(gam_c)
    Su = np.array(Su)
    # Convert all values to metric
    if units == 1:
        t = t * 0.3048  # Convert from ft to m
        Vs = Vs * 0.3048  # ft/s to m/s
        unit_wt = unit_wt * 9.81 / 62.4  # pcf to kN/m^3

    if Profile_name == '':
        Profile_name = os.path.basename(csvprofile)[:-4]

    # Create the SoilProfile object.
    sprofile = SoilProfile(Profile_name, unit_wt, Vs, t,
                          ko=ko, phi=phi, PI=PI, N160=N160,
                          layerno=layerno, layername=layername,
                          gwt=gwt, loc_mot=loc_mot,
                          outcrop_mot=outcrop_mot,
                           filename=csvprofile, OCR=OCR,
                           DS_soil=DS_soil, DS_N=DS_N, DS_freq = DS_freq,MRD = MRD,Cu = Cu,D50 = D50,gam_c = gam_c,Su = Su)
    return sprofile


def inputread(profile):
    '''
    Brings in the layer info by using a .csv input file.

    Look for the template 'TestingProfile2.csv'.   This template also
    specifies some output options. It was modelled after the SHAKE input
    option, but I decided that I didn't like them combined.  This function has
    been left for legacy reasons, but I suggest using the input_soil_profile
    function.

    The profile input file is the only input.

    The function returns the following to the SoilProfile object:
        layerno, layername, unit_wt, t, ko, phi, Vs, PI, N160, AccTHo6,
            AccTHotype6, TauTH7a, GamTH7b, RespSpec9,
            gwt, g, scalefactor, loc_mot, outcrop_mot, strainratio,
                  and Profile_name.  The AccTHo6, AccTHotype6, TauTH7a,
                  GamTH7b, and RespSpec9 variables have not been implemented
                  in the rest of the EquivLin code.

    Most of these are numpy arrays.
    The template allows you to specify the units.  This function will convert
    everything to SI units.
    '''
    # Initialize the arrays
    layerno = []
    layername = []
    unit_wt = []
    t = []
    ko = []
    phi = []
    Vs = []
    PI = []
    N160 = []
    AccTHo6 = []
    AccTHoutcrop6 = []
    AccTHotype6 = []
    TauTH7a = []
    GamTH7b = []
    RespSpec9 = []
    OCR = []
    layers = []
    openlayers = 0
    number = 0

    # Open and read the input file
    filein = open(profile, "rU")
    reader = csv.reader(filein)
    for line in reader:
        number+=1
#        print(number)
        if openlayers == 1:
            layers.append(line)
            layerno.append(int(line[0]))
            layername.append(line[1])
            unit_wt.append(float(line[2]))
            ko.append(float(line[4]))
            Vs.append(float(line[6]))

            try:
                t.append(float(line[3]))
            except:
                openlayers = 0
                continue
            phi.append(float(line[5]))
            PI.append(float(line[7]))
            N160.append(float(line[8]))
            AccTHo6.append(int(line[9]))
            AccTHoutcrop6.append(int(line[10]))
            AccTHotype6.append(int(line[11]))
            TauTH7a.append(int(line[12]))
            GamTH7b.append(int(line[13]))
            RespSpec9.append(int(line[14]))
            try:
                OCR.append(float(line[15]))
            except (IndexError,ValueError):
                pass

        if line[0] == 'Layer number':
            openlayers = 1
        elif line[0] == 'Name of Profile:':
            Profile_name = line[1]
            units = bool(line[4] == '1')
        elif line[0] == 'Scale Factor for Motion':
            scalefactor = float(line[1])
        elif line[0] == 'GWT':
            gwtlayer = int(line[1])
        elif line[0] == 'Location of Motion:':
            loc_mot = int(line[1])
        elif line[0] == 'Outcrop Motion?':
            outcrop_mot = int(line[1]) # 0 for outcrop, 1 for within
        #elif line[0] == 'Save the strain-compatible soil properties?':
        #elif line[0] == 'No. of Iterations:': iters=int(line[1])
        elif line[0] == 'Strain Ratio (for iterating)':
            try:
                strainratio = float(line[1])
            except ValueError:
                strainratio = str(line[1])
        elif line[0] == '':
            openlayers = 0
    filein.close()
    del line, openlayers

    # Convert to numpy arrays
    layerno = np.array(layerno)
    unit_wt = np.array(unit_wt)
    ko = np.array(ko)
    Vs = np.array(Vs)
    TauTH7a = np.array(TauTH7a)
    GamTH7b = np.array(GamTH7b)
    t = np.array(t)
    phi = np.array(phi)
    PI = np.array(PI)
    N160 = np.array(N160)
    AccTHo6 = np.array(AccTHo6)
    AccTHoutcrop6 = np.array(AccTHoutcrop6)
    AccTHotype6 = np.array(AccTHotype6)
    RespSpec9 = np.array(RespSpec9)
    OCR = np.array(OCR)

    # Convert all values to metric
    g = 9.81  # m/s^2
    if units == 1:
        t = t * 0.3048  # Convert from ft to m
        Vs = Vs * 0.3048  # ft/s to m/s
        unit_wt = unit_wt * 9.81 / 62.4  # pcf to kN/m^3
    gwt = sum(t[:gwtlayer - 1]) + 0.0  # Now a depth, not the top of the layer
    sprofile = SoilProfile(Profile_name, unit_wt, Vs, t,
                          ko=ko, phi=phi, PI=PI, N160=N160,
                          AccTHo6=AccTHo6, AccTHotype6=AccTHotype6,
                          AccTHoutcrop6=AccTHoutcrop6,
                          layerno=layerno, layername=layername,
                          TauTH7a=TauTH7a, GamTH7b=GamTH7b, RespSpec9=RespSpec9,
                          gwt=gwt, g=g, scalefactor=scalefactor,
                          loc_mot=loc_mot, outcrop_mot=outcrop_mot,
                          strainratio=strainratio, filename=profile, OCR=OCR)
    return sprofile


def getInMotion(motion_name, scalefactor=1.0,Butrwrth=False, b_order=4,
            cutofffreq=100, g=9.81, directout = True,**kwargs):
    '''
    Reads in the input motion from a PEER-formatted file.

    Also converts it
    to the frequency domain and applies scaling.  Also, applies a Butterworth
    lowpass filter, if specified.

    formatin = 'PEER' [default], 'PEER-OLD', or 'PEER-scaled'

    This is almost the same function as the GroundMotion.getInMotion function.
    I thought it would be useful to have it outside the object.
    '''
    formatin = kwargs.get('formatin','PEER')
    # Import the ground motion
    count = 0
    eqmot = []
    f = open(motion_name, 'rb')
    if formatin.lower() == 'peer-scaled':
       for line in f:
            count += 1
            if count == 7:
                NPTS = int(line.split()[1][:-1])
                dt = float(line.split()[3][:-4])
            if count > 7:
                for num in line.split():
                    eqmot.append(float(num))
    elif formatin.lower() == 'peer-old':
        for line in f:
            count += 1
            if count == 4:
                line = str(line)
                NPTS = int(line.split()[1][:-1])
                dt = float(line.split()[3])
            if count > 4:
                for num in line.split():
                    eqmot.append(float(num))
    #Added in by Tat as part of Tyler's request
    elif formatin.lower() == "custom":      #This gm file is basically produced using Excel
        NPTS = 0
        time = []
        for line in f:
            if NPTS == 0:
                time.append(float(line.split()[0]))
            if NPTS == 1:
                time.append(float(line.split()[0]))   
            eqmot.append(float(line.split()[1]))
            NPTS += 1
        dt = time[1]-time[0]
    else:
        for line in f:
            count += 1
            if count == 4:
                NPTS = int(line.split()[0])
                dt = float(line.split()[1])
            if count > 4:
                for num in line.split():
                    eqmot.append(float(num))
    motion = np.array(eqmot)
    if max(motion.shape) != NPTS :
        print('NPTS discrepancy! {}  {}'.format(max(motion.shape), NPTS))
        print(motion_name)
        NPTS = max(motion.shape)

    motion = motion.reshape([NPTS])
    f.close()
    del eqmot, count
    try:
        del num
    except UnboundLocalError:
        del time

    NFFT = 2 ** nextpow2(NPTS)
    f = np.linspace(0, 1/dt/2, NFFT//2+1)
    w = f[:NFFT//2+1] * 2 * np.pi
    if Butrwrth == True:
        w_cutoff=cutofffreq * 2 * np.pi
        butter = (1 / (1 + (np.abs(w) / w_cutoff)) ** (2 * b_order)) ** 0.5
    else:
        butter =  np.ones(w.shape)
    FA = np.fft.rfft(motion * scalefactor * g, NFFT) * butter
    if directout == True:
        return motion, NFFT, FA, f, w, NPTS, dt
        pass
    else:
        pass


def csvout(filename,data,columnheaders):
    '''
    This function creates a pandas dataframe, then exports it to a csv file.

    Positional Arguments:
        filename        name of output file without the .csv extension.

        data            an array of data in columns

        columnheaders   names of the columns in data

    Nothing is returned, but an output file is created.
    Requires the pandas library.
    '''
    import pandas as pd
#    pdb.set_trace()
    df = pd.DataFrame(data=data, columns=columnheaders)
    df.to_csv(filename+'.csv', index=False)


def calcAriasRatio(motion,dt, **kwargs):
    '''
    Calculates the following:
    a_rms / (max(abs(motion)) )

    It should be unitless.

    Although it uses the rms acceleration, any type of motion can be used.

    Positional Arguments:
        motion      can be a 2D array of motions (each motion is a row)
                        or a 1D array of motions

        dt          the time step of the motion(s) in seconds, can be an array

    Keyword Arguments:
        returnall   If True, returns Ar and PeakTimeRatio, else, returns only
                        Ar.  Default is False.
        duration    '5to95' [Default] or 'total'.  If '5to95' uses from 5% to
                        95% of the motion squared (Husid?) plot.
                        This keyword is passed directly to the calcRmsAccel
                        function.
    '''
    returnall = kwargs.get('returnall',False)

    if 1 in motion.shape or len(motion.shape) == 1:
        NumMotions = 1
    else:
        NumMotions = np.min(motion.shape)

    if np.isscalar(dt) & NumMotions!=1:
        dt = np.ones([NumMotions], dtype=float) * dt
    if NumMotions == 1:
        if np.isscalar(dt):
            pass
        else:
            pdb.set_trace()
            dt = dt[0]

    Ar = np.zeros([NumMotions,], dtype=float)
    PTR = np.zeros([NumMotions,], dtype=float)
    if returnall:
        if NumMotions == 1:
            RMS, PTR = calcRmsAcc(motion, dt, **kwargs)
            Ar = (RMS / np.max(np.abs(motion))) # Assume data is in rows
        else:
            for i in range(NumMotions):
#                pdb.set_trace()

                RMS, PTR[i] = calcRmsAcc(motion[i,:], dt, **kwargs)
                Ar[i] = (RMS / np.max(np.abs(motion[i,:]))) # Assume data is in rows
        return Ar, PTR
    else:
        if NumMotions == 1:
            Ar = (calcRmsAcc(motion, dt, **kwargs) /
                        np.max(np.abs(motion))) # Assume data is in rows
        else:
            for i in range(NumMotions):
                pdb.set_trace()
                Ar[i] = (calcRmsAcc(motion[i,:], dt[i], **kwargs)[0] /
                        np.max(np.abs(motion[i,:]))) # Assume data is in rows
        return Ar


def calcRmsAcc(motion, dt, **kwargs):
    '''
    Calculates the rms acceleration.

    Positional Arguments:
        motion      a 1D array of a motion, for now.

        dt          the time step of the motion(s) in seconds, a scalar

    Keyword Arguments:
        duration         '5to95' [default]  Uses the duration from 5-95% of
                            the area under the motion squared
                        'total'        Uses the entire time of the input motion

        returnall   If True, returns Ar and PeakTimeRatio, else, returns only
                        Ar.  Default is False.

    '''
    returnall = kwargs.get('returnall',False)
    motion = np.ravel(motion)
    durationtype = kwargs.get('duration','5to95')
    NPTS = len(motion)
    if np.isscalar(dt):
        pass
    else:
        dt = dt[0]
    time = np.linspace(0,(NPTS-1) * dt,NPTS)
    asqr = np.cumsum(motion ** 2) * dt
    indPeak = np.where(np.abs(motion) == np.max(np.abs(motion)))

    if durationtype == '5to95':
        int1 = interp1d( asqr, time, kind='linear')
        T = int1(0.95 * asqr[-1]) - int1(0.05 * asqr[-1])
        PeakTimeRatio = (time[indPeak] - int1(0.05 * asqr[-1])) / T
    else:
        T = time[-1]
        PeakTimeRatio = time[indPeak] / T

    if returnall:
        return np.sqrt(asqr[-1] / T), PeakTimeRatio
    else:
        return np.sqrt(asqr[-1] / T)


def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    import scipy.stats
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2



#def plotMany(NoOfPlots, xarray, yarray, specialsyms='', plottype='semilogx'):
#    ''' This is a quick way to plot several different plots at once.'''
#    if len(xarray.shape) == 1:
#        xarray=np.ones(yarray.shape)*xarray[np.newaxis,:]
#    for i in range(NoOfPlots):
#        plt.figure(i)
#        if specialsyms == '':
#            if plottype == 'semilogx':
#                plt.semilogx(xarray[i],yarray[i,:])
#            elif plottype == 'semilogy':
#                plt.semilogy(xarray[i],yarray[i,:])
#            else:
#                plt.plot(xarray[i],yarray[i,:])
#        else:
#            if plottype == 'semilogx':
#                plt.semilogx(xarray[i],yarray[i,:],specialsyms)
#            elif plottype == 'semilogy':
#                plt.semilogy(xarray[i],yarray[i,:],specialsyms)
#            else:self.MotObj.g
#                plt.plot(xarray[i],yarray[i,:],specialsyms)



#def shakeFAS(outfile1):
#    f = open(outfile1,'rb')
#    FASShake = []
#    switch = 0
#
#    for line in f:
#        if line.split() == []:
#            switch = 0
#        if switch == 1:
#            FASShake.append(line.split())
#        if line.split() == ['FREQ','FOURIER', 'AMPLITUDES']:
#            switch = 1
##            pdb.set_trace()
#    FASShake = np.array(FASShake)
#    FASShake = FASShake[:,1:]
#    FASShake = np.array(FASShake, dtype=float)
##    FASShake[:,3] = FASShake[:,0]/FASShake[:,1]
##    pdb.set_trace()
#    f.close()
#    return FASShake



#NGA_no_623_A-H05270.AT2
#NGA_no_41_ORR021.AT2
#NGA_no_1089_5081-360
#NGA_no_828_PET000.AT2
# p = "PPHS_Profile.csv"
# g = 'NGA_no_0813_YBI090.AT2'
p = "BKAN.csv"
g = "m6.0r3.4z3.0n159_84_c.AT2"
def run(profile=p, motion_name=g,
                     Degrad='DS', iters=20, error=0.01, plotting=False):
    '''
    A function that shows how to use the preceding objects.
    '''
     # Create a soil profile object.
#    sprofile = inputread(profile)  # This is the old way
    sprofile = input_soil_profile(profile)  # Better way.
#    sprofile.calcDepth(np.array([10.,25.,50.]),stresstype='mean')

    # Create a ground motion object.
    gmotion = GroundMotion(filename=motion_name,Butrwrth=True)

    # Combine the ground motion and soil profile objects, and go ahead and
    # run all the main calcs.
    result65 = Dynamics(sprofile,gmotion,Error=error,modtype='FreqInd',
                        DegradCurves=Degrad,
                        verbose=True, run_all=True,
                        strainratio='dissen') # I'm letting the script calc my dissipated energy.

    # Let's get the shear stress and strain time histories for all mid-layer depths
    taut, gamt = result65._calcTauGam(MaxOnly=False, domain='time',
                                          returns = True)

    plt.plot(gamt[0,:],taut[0,:])
    plt.xlabel('Strain of surface layer')
    plt.ylabel('Stress of surface layer')
    plt.title('Stress-Strain Hysteresis Loops')
    plt.show()

    acc = result65.calcAcc(np.sum(sprofile.t)+0.000001,domain='time',MotType='outcrop')
    plt.plot(gmotion.time,gmotion.motion,'r', label='Input')
    plt.plot(gmotion.time,acc,'b', label='Surface')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration (g)')
    plt.title('Surface and Input Acceleration Time Histories')
    plt.legend(loc='best')
#    plt.show()

    # Let's calculate the number of equivalent cycles.
    result65.calcNeq()
    result2 = Dynamics(sprofile,gmotion,Error=error,modtype='FreqInd',
                       DegradCurves=Degrad,
                       verbose=True,
                       strainratio=0.65)

    N_eq, tau_avgNeq, ratio, DE = result65.calcNeq(returns=True)

    N_eq2, tau_avgNeq2, ratio2, DE2 = result2.calcNeq(set_ratio=0.65,returns=True)
    N_eq3, tau_avgNeq3, ratio3, DE3 = result2.calcNeq(returns=True)

    # Some figures:
    fig1 = result65.plot('maxacc', syms1='-r')
    fig1 = result2.plot('maxacc', syms1='-b')

    # fig2 = result65.plot(respspec=[0.0,0.05,'within'],save=False,extension='.pdf')
    # result65.plot(respspec=[0.0,0.1,'within'],save=False,extension='.pdf',
    #                        fig=fig2 )
    # plt.figure(fig2.number)
    # plt.legend()
    # plt.show()
    # pdb.set_trace()




    result2.calcNeq()
    plt.figure('Neq')
    plt.plot(result2.N_eq,sprofile.t_mid,'-r')
    plt.plot(result65.N_eq,sprofile.t_mid,'-b')
    plt.xlabel('Number of Equiv. Cycles')
    plt.ylabel('Depth (m)')

    plt.figure('Acc')
    acc65=result65.calcAcc(MaxOnly=True)
    depths = sprofile.t_mid
    acc2 = result2.calcAcc(MaxOnly=True)
    plt.plot(acc65,depths,'-b',label='FreqInd')
    plt.plot(acc2,depths,'-r',label='FreqDep')
    plt.legend()

    FAS65 = result65.calcAcc(domain='freq',depths=[3.0])
    FAS2 = result65.calcAcc(domain='freq',depths=[3.0])
    plt.figure('FAS')
    plt.semilogx(gmotion.f,np.abs(FAS65),'-b')
    plt.semilogx(gmotion.f,np.abs(FAS2), '-r')
    plt.legend(loc='upper left')

    plt.show()

    print(result65.N_eq)

    print(result2.N_eq)
    result2.StressRedCoef()
    result65.StressRedCoef()
    plt.plot(result2.rd,sprofile.t_mid,'-r')
    plt.plot(result65.rd,sprofile.t_mid,'-b')
    plt.show()
    print(result2.D)
    print(result2.G/result2.Gmax)
    tf = result65.TF(1,0)

    plt.figure('norm')
    plt.plot(result65.gam,result65.Gratio[0,:],'-+b')
    plt.figure('log')
    plt.semilogx(result65.gam,result65.Gratio[0,:],'-+b')
    plt.show()
    # pdb.set_trace()


#    SurAcc65 = result65.calcDisp(0.1)
#    SurAccde = resultde.calcVel(MaxOnly=True)

    Sa,w = result65.RespSpec()

    plt.semilogx(w,Sa[0],'b')
    plt.semilogx(w,Sa[1],'r')



def run2():
    '''
    Asks the user via wxdialog boxes for profiles and input motions.
    This is only a start, needs to be expanded.
    '''
#    import wxfiledialog as fd

    soilfile = fd.fileopen(message = 'Select a soil profile....',
             filters = 'Comma-Separated-Values (*.csv)|*.csv',
            defaultDir='/media/Storage/Documents/Python/Modules')
    motionfile = fd.fileopen(message = 'Select a PEER motion....',
             filters = 'PEER Ground Motion (*.at2)|*.AT2||*.at2',
            defaultDir='/media/Storage/Documents/Python/StrataComp')

    sprofile = input_soil_profile(soilfile[0])
    Vs30, H = sprofile.calcVs30(returns=True)
    Tsoil = 4. * H / Vs30
    print('{} {}'.format(Tsoil, H))
    gmotion = GroundMotion(filename=motionfile[0],Butrwrth=False)

    result = Dynamics(sprofile,gmotion,modtype='FreqInd',
                        verbose=True)
    Sa,w = result.RespSpec(depths=[0.])
    Sa2,w2 = result.RespSpec(depths=[0.],w=[2*np.pi*Vs30/(4*H)])
    print('{}  {}'.format(Sa2, w2))
#    pdb.set_trace()

#    csvout('test',np.column_stack((2*np.pi/w,Sa)),['Period','SpectralAcceleration'])
    plt.semilogx(2*np.pi/w,Sa)
#    sprofile.plot()
#    result.plot('all')
    plt.show()
    pdb.set_trace()


def run3():
    import wxfiledialog as fd


    motionlist = fd.fileopen(message = 'Select a PEER motion....',
             filters = 'PEER Ground Motion (*.at2)|*.AT2||*.at2',
            defaultDir='/media/Storage/Documents/00-VtResearch/My Publications/10NCEE/Comparisons/Simple/Motions/')
#    motionlist = ['NGA_no_0813_YBI090.AT2','NGA_no_623_A-H05270.AT2']
    Ar = []
    Ar2 = []
    pf = []
    Neq = []
    for motion in motionlist:
        sprofile = SoilProfile('ThinSoil',
                               unit_wt = np.array([18,22], dtype=float),
                               Vs = np.array([250,650], dtype=float),
                                t = np.array([1.], dtype=float),
                                ko = np.array([0.4]) )
        gmotion = GroundMotion(name = motion)
    #    gmotion = GroundMotion(name = 'NGA_no_623_A-H05270.AT2')
        # Get predomninant Freq:
#        [freq for (fas,freq) in sorted(zip(gmotion.f,np.abs(gmotion.FA)))]
#        fas,freq = sorted(zip(gmotion.f,np.abs(gmotion.FA)))


        result = Dynamics(sprofile, gmotion, verbose=False, run_all=True)

#        print('N_eq: {}'.format(result.N_eq))
#        print('N_eq: {}'.format(result.calcNeq(depths=[1.001],returns=True)[0]))
#        print('Ar: {}'.format(gmotion.Ar))
#        print('Predom Freq: {}'.format(gmotion.predfreq))
#        plt.semilogx(gmotion.f,np.abs(gmotion.FA))
#        plt.show()
#        print('Ar2: {}'.format(gmotion.Ar2))
#        pdb.set_trace()
        Ar.append(gmotion.Ar)
        Ar2.append(gmotion.Ar2)
        pf.append(gmotion.predfreq)
        Neq.append(result.calcNeq(depths=[1.001],returns=True)[0])

    Ar = np.array(Ar)
    Ar2 = np.array(Ar2)
    pf = np.array(pf)
    Neq = np.array(Neq).ravel()
    pdb.set_trace()


def run4(profile='Profile020.csv', motion_name='NGA_no_12_PEL090.AT2',
                     Degrad='IZ', iters=60, error=0.03, plotting=False):
    '''
    A function that shows how to use the preceding objects.
    '''
     # Create a soil profile object.
#    sprofile = inputread(profile)  # This is the old way
    sprofile = input_soil_profile(profile)  # Better way.
#    sprofile.calcDepth(np.array([10.,25.,50.]),stresstype='mean')
    pdb.set_trace()
    # Create a ground motion object.
    gmotion = GroundMotion(filename=motion_name,Butrwrth=False,
            formatin='PEER')

    # Combine the ground motion and soil profile objects, and go ahead and
    # run all the main calcs.
    result65 = Dynamics(sprofile,gmotion,Error=error,modtype='FreqInd',
                        DegradCurves=Degrad,
                        iters=iters,
                        verbose=True, run_all=True,
                        strainratio=0.65)#'dissen') # I'm letting the script calc my dissipated energy.
    tau = result65.calcStress(2, domain='time')
    gam = result65.calcStrain(2, domain='time')
    time = gmotion.time

    de = result65.calcDissEn(tau.reshape(1, len(tau)), gam.reshape(1, len(tau)), FinalValOnly=False).T
    ind = np.where(np.sign(tau[1:]) != np.sign(tau[:-1]))
    plt.plot(time[:-1], de)
    plt.plot(time[:-1][ind], de[ind])
    plt.show()
    pdb.set_trace()

# if __name__ == '__main__':
    # run()
#    run2()
#    run3()
#    run4()
