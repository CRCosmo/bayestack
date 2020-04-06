#!/usr/bin/env python

"""
This is plot_counts.py
Jonathan Zwart
May 2014

Usage:

./plot_counts.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
if os.getenv('PBS_O_HOST') not in ['baltasar']:
    from matplotlib import pyplot as plt
import pylab
from profile_support import profile
import pymultinest
import lumfuncUtils
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats,get_changed_multiplot
from countUtils import calculateDnByDs,medianArray,powerLawFuncQuadS
from lumfuncUtils import get_Lbins,get_sbins,get_dl, get_Vmax,get_dsdl,get_z
from plat import open_anytype, calc_Vmax, calcu_zeff,lumfuncs
from matplotlib.ticker import AutoMinorLocator
from cycler import cycler
#import sourcecounts1

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

@profile
def main():

    """
    """

    print 'Settings file is %s' % setf

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
    recon_expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
    #print expt.data
    #print expt.bins

    # Get MAP parameters
    print '\n'*5
    print len(expt.parameters)
    ncols = len(expt.parameters)
    summf=os.path.join(outdir,'1-summary.txt')
    summary=numpy.genfromtxt(summf)[-1,:]
    drawmap=summary[-(ncols+2):-2] 
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawmap=ana.get_best_fit()['parameters']
    
    if dataset=='sdss':
        print datafiles, len(datafiles[0])
        dataf = datafiles[0]
        if len(dataf)==20:
            num = dataf[-5]
            print num
        else:
            num = dataf[-6:-4]
        d=numpy.genfromtxt('%s/%s'%(outdir,datafiles[0][4:]))
        print datafiles[0][4:]
    elif dataset=='cosmos':
        print datafiles, len(datafiles[0])
        dataf = datafiles[0]
        if len(dataf)==24:
            num = dataf[-5]
            print num
        else:
            num = dataf[-6:-4]
        d=numpy.genfromtxt('%s/%s'%(outdir,datafiles[0][8:]))
        print datafiles[0][8:]

    elif 'sim' in dataset:
        d=numpy.genfromtxt(os.path.join(outdir,'sim.txt'))
    
    #bins2 = numpy.arange(21.4,29.2,0.4)
    bins2 = numpy.arange(18.0,29.2,0.4)
    bins3 = numpy.arange(18.,29.2,0.01)
    print len(expt.parameters)
    print 'this is 1st num ', num
    #sys.exit()
    params = drawmap# [0]
    #for i in drawmap:
    #	params.append(i)
    print params
    print expt.parameters
    #novak 2017
    L_n   =[21.77,22.15,22.46,22.77,23.09,23.34]
    rho_n =[-2.85,-2.88,-3.12,-3.55,-4.05,-4.63]
    L_ner_u =[0.23,0.18,0.19, 0.2, 0.21, 0.28]
    L_ner_d =[1.1, 0.15, 0.14, 0.12, 0.12, 0.048]
    rho_ner=[0.09, 0.03,0.0355,0.059,0.1,0.234]

     #novak 2018
    L_n2   =[[21.77,22.24,22.68,23.16,23.69,24.34, 24.74,25.56],\
             [22.30,22.61,22.96,23.38,23.80,24.10, 24.55,25.14],\
             [22.61,22.86,23.14,23.45,23.82,24.14, 24.40,24.71],\
             [22.85,23.16,23.69,24.24,24.80,25.31, 25.96,26.69],\
             [23.10,23.38,23.86,24.36,24.86,25.35, 25.94,26.36],\
             [23.32,23.57,23.94,24.32,24.67,25.06, 25.47,25.96],\
             [23.54,23.75,24.13,24.45,24.90,25.27, 25.74,26.10],\
             [23.73,23.99,24.57,25.10,25.68,26.18, 26.83,27.51],\
             [24.01,24.26,24.76,25.26,25.91,26.19, 27.45],\
             [24.30,24.56,24.80,25.13,25.37,25.80, 25.97,26.49]]
             
 
    L_ner_u2  =[[0.23,0.27, 0.34, 0.37, 0.35, 0.21, 0.31, 0.03],\
                [0.11,0.20, 0.26, 0.25, 0.24, 0.35, 0.31, 0.15],\
                [0.08,0.15, 0.20, 0.22, 0.17, 0.18, 0.24, 0.28],\
                [0.072,0.33,0.38, 0.40, 0.48, 0.42, 0.41, 0.28],\
                [0.080,0.30,0.32, 0.32, 0.32, 0.33, 0.24, 0.34],\
                [0.068,0.22,0.24, 0.25, 0.29, 0.30, 0.28, 0.21],\
                [0.067,0.22,0.21, 0.26, 0.18, 0.17,0.075, 0.11],\
                [0.093,0.36,0.31, 0.30, 0.25, 0.28, 0.16,0.028],\
                [0.076,0.34,0.36, 0.38, 0.24, 0.47, 0.27],\
                [0.097,0.14,0.20, 0.16, 0.23, 0.092,0.22,0.026]]
                
    
    L_ner_d2  =[[1.0, 0.24, 0.17, 0.14, 0.16, 0.30, 0.20, 0.50],\
               [0.29, 0.21, 0.15, 0.16, 0.17, 0.06, 0.099,0.29],\
               [0.24, 0.17, 0.12, 0.10, 0.15, 0.15, 0.08,0.074],\
               [0.16, 0.24, 0.19, 0.17, 0.15, 0.09, 0.16, 0.33],\
               [0.30, 0.19, 0.18, 0.18, 0.18, 0.17, 0.26, 0.18],\
               [0.14, 0.18, 0.15, 0.14, 0.10, 0.095,0.11, 0.21],\
               [0.19, 0.14, 0.15, 0.11, 0.19, 0.19, 0.29, 0.28],\
               [0.40, 0.17, 0.21, 0.23, 0.27, 0.25, 0.37, 0.53],\
               [0.20, 0.17, 0.16, 0.14, 0.27, 0.042, 0.27],\
               [0.22, 0.16, 0.10, 0.14,0.072, 0.21, 0.08, 0.30]]

    rho_n2 =[[-2.84,-2.90,-3.34,-4.00,-4.92,-4.92, -5.22,-5.07],\
             [-2.95,-3.17,-3.46,-4.24,-4.76,-5.41, -5.23,-5.44],\
             [-2.95,-3.11,-3.44,-3.81,-4.37,-4.48, -4.90,-5.14],\
             [-2.99,-3.24,-3.89,-4.54,-5.27,-5.33, -5.69,-5.89],\
             [-3.32,-3.51,-4.06,-4.74,-5.28,-5.43, -6.08,-5.70],\
             [-3.29,-3.54,-4.02,-4.45,-5.11,-5.56, -5.37,-6.07],\
             [-3.43,-3.61,-4.19,-4.56,-5.09,-5.34, -5.70,-5.75],\
             [-3.88,-3.98,-4.64,-5.33,-5.73,-6.27, -6.14,-6.77],\
             [-3.85,-4.28,-5.05,-5.60,-5.68,-6.44, -6.59],\
             [-4.43,-4.91,-5.46,-5.58,-5.89,-6.23, -6.26,-6.91]]
       
    rho_ner_u2=[[0.08,0.024,0.038,0.083,0.28,0.28 , 0.45, 0.34],\
                [0.048,0.027,0.036,0.089,0.17,0.45, 0.34, 0.45],\
                [0.062,0.024,0.032,0.048,0.096, 0.11,0.18,0.25],\
                [0.042,0.018,0.035,0.073,0.180, 0.20,0.34,0.45],\
                [0.045,0.019,0.033,0.072,0.15,  0.17,0.45,0.25],\
                [0.040,0.022,0.033,0.055,0.13,  0.22,0.17,0.45],\
                [0.061,0.020,0.035,0.055,0.099, 0.14,0.22,0.25],\
                [0.044,0.022,0.044,0.11, 0.16, 0.34, 0.28,0.76],\
                [0.048,0.025,0.058,0.13, 0.24, 0.34, 0.45],\
                [0.084,0.058,0.11 ,0.24, 0.20, 0.28, 0.28,0.76]]
                
    rho_ner_d2=[[0.07,0.023,0.035,0.070,0.25,0.25 , 0.37, 0.30],\
                [0.043,0.026,0.033,0.074,0.16,0.37, 0.30, 0.37],\
                [0.054,0.023,0.030,0.043,0.079,0.086,0.17,0.22],\
                [0.038,0.017,0.032,0.062,0.17,0.19, 0.30, 0.37],\
                [0.041,0.019,0.031,0.062,0.11,0.16, 0.37, 0.22],\
                [0.037,0.021,0.031,0.049,0.098,0.20,0.16, 0.37],\
                [0.054,0.019,0.032,0.049,0.081,0.11,0.20, 0.22],\
                [0.040,0.021,0.040,0.085,0.15, 0.30, 0.25,0.52],\
                [0.043,0.024,0.051,0.10, 0.15, 0.30, 0.37],\
                [0.070,0.051,0.087,0.16,0.19, 0.25, 0.25, 0.52]]

    
    fig = plt.figure()
    pre_chain='01'
    #pre_chain='01'
    chains=['12a','12b','12c','12d','12e','12f','12g','12h','12i']#chains_191062
    chains=['22a','22b','22c','22d','22e','22f','22g','22h','22i']
    chains=['32a','32b','32c','32d','32e','32f','32g','32h','32i']
    chains=['42a','42b','42c','42d','42e','42f','42g','42h','42i']
    chains=['62a_2','62b_2','62c_2','62d_2','62e_2','62f_2','62g_2','62h_2','62i_2'] #chains_191062
    chains=['01a_3','01b_3','01c_3','01d_3','01e','01f_3','01g_3','01h_3','01i_3']#chains_191101
    chains=['01a_4','01b_4','01c_4','01d_4','01e_4','01f_4','01g_4','01h_4','01i_4_2']#chains_191101 7x7 grid
    chains=['02a_4','02b_4','02c_4','02d_4','02e_4','02f_4','02g_5','02h_4','02i_4']#chains_191101 7x7 grid lognorm
    #chains=['01a_5','01b_5','01c_5','01d_5','01e_5','01f_5','01g_5','01h_5','01i_5']#chains_191101 3x3 grid scaled
    #chains=['01a_6','01b_6_3','01c_6','01d_6_2','01e_6_2','01f_6','01g_6_2','01h_6_2','01i_6_2']#chains_191101 11x11 grid
    #chains=['01a','01b','01c','01d','01e','01f','01g','01h','01i'] #chains_1912 #includeing maskregion
    chains=['01a_2_1','01b_2_1','01c_2_1','01d_2','01e_2','01f_2','01g_2','01h_2','01i_2'] #1912 exclding peterinclde match
    #chains=['01a_2','01b_2','01c_2_1','01d_2','01e_2','01f_2','01g_2','01h_2','01i_2'] #chains_1912 excluding peter
    #chains=['01a_3','01b_3','01c_3','01d_3','01e_3','01f_3_1','01g','01h','01i'] #chains_1912 including peter
    #chains=['01a_4_1','01b_4_1','01c_4','01d_4','01e_4','01f_4_1','01g_4','01h_4','01i_4'] #chains_1912 no mass selection
    chains=['02a','02b_2','02c','02d','02e','02f','02g','02h','02i'] #lgnorm #chains_1912 excluding peter
    chains=['01a_1','01b','01c','01d','01e','01f_1','01g','01h','01i'] #chains_2001 DR4
    chains=['01a_4','01b_2','01c_2','01d_2','01e_2','01f_2','01g','01h_2','01i_2'] #chains_2001 DR4 remove phi2
    chains=['01a_8','01b_8','01c_8','01d_8','01e_8','01f_8','01g_7','01h_8','01i_8','01j_7'] #chains_2001 DR4 dpl_dpl
    chains2=['02a_z_6x','02b_z_6x','02c_z_6x','02d_z_6x','02e_z_6x','02f_z_6x','02g_z_6x','02h_z_6x','02i_z_6x','02j_z_6x']
    chains3='02z_31_2'
    #z = [0, 0.5, 1 , 1.5, 2, 2.3, 2.6, 3, 3.5, 4] #old
    #z = [ 0.7, 1 , 1.35, 1.7, 2, 2.3, 2.6, 3, 3.5, 4]# older
    print expt.parameters
    n,m=3,3
        
    #print xrecon
    #print yrecon
    #sys.exit()
    
    z_new=False
    Replot=False
    z_nov = [0.31, 0.5, 0.69, 0.9, 1.16, 1.44, 1.81, 2.18, 2.81, 3.71, 4.83]
    if z_new:
        z =[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,1.9,2.2,2.5, 2.8,3.2,3.6,4.0]
        chains=['11a','11b','11c_3','11d_3','11e_3','11f_3','11g_3','11h','11i','11j', '11k', '11l' ,'11m', '11n']#chains_201111
    else:
        z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]#oldest ;)
    
    ana3=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_2002%s'%chains3,outstem))
    drawmap3=ana3.get_best_fit()['parameters']
    params3 = drawmap3
    
    parameters2= ['noise', 'A_SF', 'A_agn', 'LMIN', 'LMAX']
    parameters3= ['noise', 'A_agn', 'A_SF', 'B_agn', 'B_SF', 'LMIN', 'LMAX']
    parameters3= ['noise', 'A_agn', 'A_SF', 'B_agn', 'B_SF', 'LMIN_1', 'LMIN_2', 'LMIN_3', 'LMIN_4', 'LMIN_5', 'LMIN_6', 'LMIN_7', 'LMIN_8', 'LMIN_9', 'LMIN', 'LMAX']
        
    for num in range(1,10):
    #num = 1.
     print num
     print 'new num ',num
     z_min,z_max, z_m = get_z(num,z_new)
     dl = get_dl(z_m)
     Vmax = get_Vmax(z_min,z_max)
     dsdl = get_dsdl(z_m,dl)

     z_m = (z_min + z_max)/2                       
     dl = get_dl(z_m)
    
     print z_min,z_max
     print expt.kind
     area = SURVEY_AREA*sqDeg2sr 
     area1 = SURVEY_AREA*sqDeg2sr
     print SURVEY_AREA
     ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_20%s%s'%(pre_chain,chains[num-1]),outstem))
     drawmap=ana.get_best_fit()['parameters']
     drawmid=ana.get_stats()['marginals']
     params = drawmap
     
     ana2=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_2002%s'%chains2[num-1],outstem))
     drawmap2=ana2.get_best_fit()['parameters']
     params2 = drawmap2
  
     
     settingf='chains_20%s%s/bayestack_settings.py'%(pre_chain,chains[num-1])
     f = open(settingf, 'r')
     mod = f.readlines()[31].split()[0]
     print mod
     if 'lognorm' in mod:
        if 'dpl' in mod:
            model='LFlognorm_dpl'
        elif 'pl_' in mod:
            model='LFpl_lognorm'
        else:
            model='LFlognorm'
     elif 'LFdpl' in mod or 'LFpl' in mod:
        if 'LFdpl_dpl' in mod:
            model='LFdpl_dpl'
        elif '_pl' in mod:
            model='LFdpl_pl'
        elif 'LFpl_' in mod:
            model='LFpl_dpl'
        elif 'LFpl' in mod:
            model='LFpl'
        else:
            model = 'LFdpl'
     
     with open('chains_20%s%s/params.tex'%(pre_chain,chains[num -1])) as f:
        parameters=[line.split(None,1)[0] for line in f]
    #model='LFdpl_pl'
     print parameters
     print model
     print 'LSIGMA' in parameters
     
     if model in ['LFsch','LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFlognorm_dpl','LFpl_lognorm']:
        Lnorm=params[parameters.index('LNORM')]
        Lstar=params[parameters.index('LSTAR')]
        Lslope=params[parameters.index('LSLOPE')]
        #Lzevol=params[parameters.index('LZEVOL')]
     if model in ['LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFlognorm_dpl']:
        Lslope2=params[parameters.index('LSLOPE2')]
     if model in ['LFlognorm','LFpl', 'LFdpl', 'LFdpl_pl','LFpl_dpl', 'LFlognorm_dpl','LFpl_lognorm','LFdpl_dpl']:
        Lnorm_2=params[parameters.index('LNORM_2')]
        Lstar_2=params[parameters.index('LSTAR_2')]
        Lslope_2=params[parameters.index('LSLOPE_2')]

     if model in ['LFdpl_dpl','LFdpl']:
        Lslope2_2=params[parameters.index('LSLOPE2_2')]
    
     if model in ['LFlognorm','LFlognorm_dpl','LFpl_lognorm']:
        Lsigma = params[parameters.index('LSIGMA')]
     
     Lmin=params[parameters.index('LMIN')]
     Lmax=params[parameters.index('LMAX')]
     
     Lmin=params[expt.parameters.index('LMIN')]
     Lmax=params[expt.parameters.index('LMAX')]

     Lmin=0.7
     SMIN  = get_sbins(numpy.power(10,Lmin),z_m,dl)*1e6
     SMAX  = get_sbins(numpy.power(10,Lmax),z_m,dl)*1e6
     sigma,fsigma,Lmin = numpy.log10(get_Lbins([SURVEY_NOISE,SURVEY_NOISE*5, Lmin],z_m,dl,'muJy')*(1.4/3)**(-.7))
     sigma2,fsigma2 = numpy.log10(get_Lbins([450,2400],z_m,dl,'muJy'))
     print z_m,dl, sigma,fsigma
     print SURVEY_NOISE,SURVEY_NOISE*5,SURVEY_AREA
     
     s=numpy.loadtxt('chains_20%s%s/recon_stats.txt'%(pre_chain,chains[num -1]))
     xrecon=s[:-1,0]; yrecon=s[:-1,1]
     yrecon_d=s[:-1,2]; yrecon_u=s[:-1,3]
     yrecon_rms = s[:-1,4]
     yrecon_avr = s[:-1,5]
     yrecon_rms_down = yrecon_rms
    
     #area = 6672*sqDeg2sr #143165.15  49996.49 59876.8861135
     sbin3 = get_sbins(10**bins2,z_m,dl)*1e6
     sbin4 = get_sbins(10**bins3,z_m,dl)*1e6
     ssd = numpy.loadtxt('chains_20%s%s/data_cos_s%s.txt'%(pre_chain,chains[num-1],num))
     ssbin1 = ssd[:,0] #2 
     ssbin  = ssd[:,2] #0
     Nbin   =  ssd[:,3] 
     zbin  = ssd[:,-1]
         
     L_s   = numpy.log10(get_Lbins(ssbin,z_m,dl,'muJy')*(1.4/3.)**(-.7))
     #print expt.binsMedian
     #sys.exit()
    
     xreco = get_Lbins(xrecon,z_m,dl,'muJy')#*(1.4/3)**(-.7)
     yreco = yrecon
    
    
    
     print SMIN,SMAX
     #sys.exit()
     lin_yrecon = numpy.log10(yreco)
     #print xrecon
     lin_xrecon = numpy.log10(xreco)
    
    
     bins2 = numpy.arange(15.,25.,0.4)
     bins  = numpy.arange(fsigma,26.2,0.4)#bin +0.2
     if not os.path.isfile('cos_data/cos_s%s_LF.el'%num) and z_new or not os.path.isfile('cos_data/cos_s%s_LF_old.el'%num) and not z_new or Replot:
      print 'doing all calculations'
      #l,z_l = open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,21), F='uJy',L=True,getz=True)
      #l2,z_l2 = open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,24), F='uJy',L=True,getz=True)
      if z_new:
        L,z_L = open_anytype('cos_data/cos_s%s.txt'%num,(2,3), F='uJy',L=True,getz=True,band='S')
        print 'this is the new z bins (please check if lumfuncUtils and cos_manifest are using the right z)'
      else:
        L,z_L = open_anytype('cos_data/cos_s%s_old.txt'%num,(2,3), F='uJy',L=True,getz=True,band='S')
        print 'this is the old z binning style. (please check if lumfuncUtils and cos_manifest are using the right z)'
      L_c,z_c=open_anytype('cosmos_counter_parts_II.txt',(3,4),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
      z_c2,S_c2,L_c2 =numpy.loadtxt('cosmos_counter_parts_II.txt',unpack=True, usecols=(3,4,5))
      L_c2=L_c2[(z_c2<z_max)&(z_c2 >z_min)]
      S_c2=S_c2[(z_c2<z_max)&(z_c2 >z_min)]
      z_c2=z_c2[(z_c2<z_max)&(z_c2 >z_min)]
    
    
    #rho_1,er1 = calc_Vmax(calcu_zeff(z_l,l,11.5e-32),l,bins,z_max,z_min,area)
    #rho_2,er2 = calc_Vmax(calcu_zeff(z_l2,l2,11.5e-32),l2,bins,z_max,z_min,area)
      rho_3,er3 = calc_Vmax(calcu_zeff(z_L,L,1.15e-31),L,bins,z_max,z_min,area)
    #rho_4,er4 = calc_Vmax(99.,L,bins2,z_max,z_min,area)
      rho_5,er5 = calc_Vmax(calcu_zeff(z_c,L_c,1.15e-31),L_c,bins,z_max,z_min,area)
      rho_6,er6 = calc_Vmax(calcu_zeff(z_c2,10**L_c2,1.15e-31),10**L_c2,bins,z_max,z_min,area)
      if z_new:
        f = open('cos_data/cos_s%s_LF.el'%num,'w')
      else:
        f = open('cos_data/cos_s%s_LF_old.el'%num,'w')
      for i in range(len(bins)):
        f.write('%5.1f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n'%(bins[i], rho_3[i],er3[i],rho_5[i],er5[i],rho_6[i],er6[i]) )
      f.close()
    
      dl_c=get_dl(z_c2)
      L_c3=numpy.log10(get_Lbins(S_c2,z_c2,numpy.array(dl_c),'muJy')*(1.4/3.)**(-.7))
     else:
        if z_new:
            bins,rho_3,er3,rho_5,er5,rho_6,er6 = numpy.loadtxt('cos_data/cos_s%s_LF.el'%num, unpack=True)
            print 'Using stored RLF for the new z (The one with more z bins). If you want to calculate the plot the set "Replot=True"'
        else:
            bins,rho_3,er3,rho_5,er5,rho_6,er6 = numpy.loadtxt('cos_data/cos_s%s_LF_old.el'%num, unpack=True)  
            print 'Using stored RLF for the is old z bins. If you want to calculate the plot the set "Replot=True"'       
     #ssd = numpy.loadtxt('%s/data_cos_s%s.txt'%(outdir,num))
     ssd = numpy.loadtxt('cos_data/cos_s%s.txt'%(num), unpack=True, usecols=(-1,2))
     ssbin2 = ssd[:,0] #2 
     ssbin1 = ssd[:,2] #0
     Nbin2  =  ssd[:,3] 
     zbin3  = ssd[:,-1]
     dl1 = get_dl(z_min)
     dl2 = get_dl(z_max)
    #dl3 = get_dl(zbin3)
    
     Ncount = sum(Nbin2[fsigma>ssbin1])    
     Llmin = numpy.log10(get_Lbins([fsigma/numpy.sqrt(Ncount)],z_min,get_dl(z_min),'muJy'))[0]
     Llmin*=(1.4/3.)**(-.7)
     print Llmin
    #sys.exit()
        
     fig.add_subplot(n,m,get_changed_multiplot(n,m,num))
    
     plt.rc('lines', linewidth=2)
     
     #plt.errorbar(bins,rho_1 ,yerr=er1,fmt='*r',label='detected extracted',markersize=11) #bin
     
     
     
     if model in ['LFdpl','LFdpl_dpl']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
     elif model in ['LFpl','LFdpl_pl']:
        faint =lumfuncUtils.powerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lnorm_2)
     elif model in ['LFpl_dpl']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),Lstar,Lslope,Lslope2,Lnorm)
     else:
        faint = lumfuncUtils.lognormpl(numpy.power(10,bins3), Lstar_2,Lslope_2,Lsigma,Lnorm_2)
        
    
     phi_dpl =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),Lstar ,Lslope,Lslope2,Lnorm)
    
    #phi_dpl
    #dpl  = 
     dm =0.4#05
     hm=0.
     #if num==1:hm=-0.25 

     L_1 = numpy.power(10,bins3)/(1 + z_nov[num-1])**(2.86-z_m*0.7)
     L_2 = numpy.power(10,bins3)/(1 + z_nov[num-1])**(2.95-z_m*0.29)
     Novak17=lumfuncUtils.lognormpl(L_2, numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
     Novak17_2=lumfuncUtils.lognormpl(numpy.power(10,bins3)/(1 + z_nov[num-1])**(3.16-z_nov[num-1]*0.32),\
        numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3/2.5))
     phi_agn_n=lumfuncUtils.doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(10**(-5.5)))
     phi_tot_n = phi_agn_n+Novak17_2
     phi_2 =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params2,paramsList=parameters2,\
        family='LFevol_logn_mat',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True) for lL in bins3]
     phi_3 =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params3,paramsList=parameters3,\
        family='LFevol_logn',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True) for lL in bins3]
     #sys.exit()
     plt.plot(bins3+hm,numpy.log10(faint)-dm,'--b' ,label='SF MAP')
     #plt.plot(numpy.log10(10**bins3),numpy.log10(phi_dpl)-dm,'--c', label='agn')
     #plt.plot(numpy.log10(10**bins3 *(1.4/3)**(-.7))+hm,numpy.log10(faint)-dm,'--b' ,label='faint-end')   
     plt.errorbar(lin_xrecon+hm,lin_yrecon-dm ,fmt='-b', label='MAP',linewidth=1.8)
     
     plt.errorbar(numpy.log10(10**bins3/(1.4/3)**(-.7)), numpy.log10(phi_3) -0.4,fmt='.g',linewidth=1.2,label='Evolution Novak LF Model fit')
     plt.errorbar(numpy.log10(10**bins3), numpy.log10(phi_2) -0.4,fmt='-.k',linewidth=4.5,label='Indvidual Novak LF Model fit')
     #plt.errorbar(lin_xrecon+hm,numpy.log10(yrecon_avr)-dm,fmt='-c', label = 'Average')
     plt.errorbar(numpy.array(L_n2[num-1]),numpy.array(rho_n2[num-1]) -0.4,xerr=[L_ner_d2[num-1],L_ner_u2[num-1]],yerr=[rho_ner_d2[num-1],rho_ner_u2[num-1]], fmt='sr',fillstyle='none',ecolor='k',markeredgewidth=1.5, label='Total RLF Novak et al 2018')
     plt.errorbar(numpy.log10(10**bins3), numpy.log10(phi_tot_n),fmt='--r',label='Total fit Novak et al 2018')
     #plt.errorbar(bins3, numpy.log10(Novak17)-0.4,fmt='--b',label='SF galaxies Novak et al 2017')
     #plt.errorbar(bins3, numpy.log10(Novak17_2),fmt='--r',label='SF galaxies Novak et al 2017')
     #plt.plot(bins3,numpy.log10(phi_agn_n),'--g', label='local agn nov')
     
     #plt.plot(L_s, numpy.log10(rho_1), '*m',label='bins LF',markersize=8)
     #plt.errorbar(bins,rho_5 ,yerr=er5,fmt='*r',label='Total',markersize=10) #bin
     #plt.errorbar(bins,rho_6 ,yerr=er6,fmt='>g',label='Total LF',markersize=10) #bin
     plt.fill_between(lin_xrecon+hm,numpy.log10((yrecon_d))-dm,numpy.log10(yrecon_u) -dm, color='k',alpha=0.2)
     plt.errorbar(bins[1:],rho_3[1:] ,yerr=er3[1:],fmt='h', fillstyle= 'none', markeredgewidth=3 ,label='extracted fluxes',markersize=12) #bin
    
     
     plt.ylim(-7.6,-1.6)
     plt.axvline(fsigma,color='k',linestyle='dashed')
     plt.axvline(sigma,color='g')
     #plt.axvline(Lmin,color='r',linewidth=2)
     print 'Is this still num ',num
     if num<2 or num<6 and num>3 or num>6 and num<9:
        plt.xticks( visible=False)
        plt.ylim(-7.6,-1.6)
     
     if num<4:
        plt.xlim(19,25.8)
        plt.text(19.2,-7.,r'$\rm{%.1f < z < %.1f}$'%(z_min,z_max),fontsize = 25 ) # \n $M_i$ < -22
     elif num<7:
        plt.xlim(20.8,27.2)
        plt.text(21.,-7.,r'$\rm{%.1f < z < %.1f}$'%(z_min,z_max),fontsize = 25 ) # \n $M_i$ < -22
     else:
        plt.xlim(21.3,27.8)   
        plt.text(21.5,-7.,r'$\rm{%.1f < z < %.1f}$'%(z_min,z_max),fontsize = 25 ) # \n $M_i$ < -22
        
    
     if num==3:
    	plt.xticks([19,20,21,22,23,24,25],fontsize=22)
     if num==6:
       plt.xticks([21,22,23,24,25,26,27],fontsize=22)
     if num == 9:
    	plt.xticks([22,23,24,25,26,27],fontsize=22)
    	plt.legend(frameon=False).draggable()
     if num==1 or num ==2 or num==3 or num==10:
         plt.yticks([-2.,-3,-4,-5, -6], fontsize=22)
         if num==3:
            plt.yticks([-2.,-3,-4,-5, -6, -7], fontsize=22)
         continue
     plt.yticks(visible=False) 
    
    fig.text(0.06,0.5,r'$\rm{log_{10}[\Phi /(Mpc^{-3}mag^{-1})]}$',fontsize=30, ha='center',va='center', rotation='vertical')
    fig.text(0.5,0.04,r'$\rm{log_{10}[L_{1.4}/(W}$ $\rm{Hz^{-1})]}$',fontsize = 30, ha='center',va='center')    
    #plt.text(23.9,-9.7,'%s'%outdir,fontsize = 16 ) # \n $M_i$ < -22
    plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    plt.tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    plt.subplots_adjust(hspace=0,wspace=0)
    #plt.xtick(fontsize=15)
    #plt.ytick(fontsize=15)
    
    #plt.ylim(-10.2,-5.8)
    #plt.ylim(-11,-5.5)
     #ax.xaxis.set_minor_locator(AutoMinorLocator())
     #ax.yaxis.set_minor_locator(AutoMinorLocator())
    #print truth['LMIN']
    
    plotf='%s/LF_recon_s1.pdf' % (outdir)
    plt.show()
    
	
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
