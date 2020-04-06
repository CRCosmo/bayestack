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
from lumfuncUtils import get_Lbins,get_sbins,get_dl, get_Vmax,get_dsdl,get_z,get_sfr_q, get_sfrd_z, get_q, sfrd_Behroozi, sfrd_Madau
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

    nov_z,nov_z_u, nov_z_l,  nov_sfrd, nov_sfrd_u, nov_sfrd_l,  Lower, lw_up, lw_lw, nov_phi_L, nov_phi_L2 = numpy.loadtxt('novak_2017.txt',usecols=(0,1,2, 3,4,5, 6,7,8,-2,-1),unpack=True)
    
    fig = plt.figure()
    pre_chain='01'
    chains=['01a_4','01b_4','01c_4','01d_4','01e_4','01f_4','01g_4','01h_4','01i_4_2', '01j_4_2']#chains_191101 7x7 grid
    chains=['02a_4','02b_4','02c_4','02d_4','02e_4','02f_4','02g_5','02h_4','02i_4','02j_4']#chains_191101 7x7 grid lognorm
    #chains=['01a_5','01b_5','01c_5','01d_5','01e_5','01f_5','01g_5','01h_5','01i_5']#chains_191101 3x3 grid scaled
    #chains=['01a_6','01b_6','01c_6','01d_6','01e_6','01f_6','01g_6','01h_6','01i_6']#chains_191101 11x11 grid
    #chains=['01a','01b','01c','01d','01e','01f','01g','01h','01i'] #chains_1912 #includeing maskregion
    chains=['01a_3','01b_3','01c_3','01d_3','01e','01f','01g','01h','01i','01j'] #chains_1912 including peter
    #chains=['01a_2_1','01b_2_1','01c_2_1','01d_2','01e_2','01f_2','01g_2','01h_2','01i_2','01j_2']
    #chains=['01a_4_1','01b_4_1','01c_4','01d_4','01e_4','01f_4_1','01g_4','01h_4','01i_4','01j_4'] #chains_1912 no mass selection
    #chains=['01a_3','01b_2','01c_2','01d_2','01e_2','01f_2','01g','01h_2','01i_2','01j_3'] #chains_2001 DR4 remove phi2
    chains=['01a_5','01b_5','01c_5','01d_5','01e_5','01f_5','01g_5','01h_5','01i_5','01j_5'] #chains_2001 DR4 remove phi2 new changes
    chains=['01a_8','01b_8','01c_8','01d_8','01e_8','01f_8','01g_7','01h_8','01i_8','01j_7'] #chains_2001 DR4 dpl_dpl
    #chains=['02a_z_9_6','02b_z_9_6','02c_z_9_6','02d_z_9_6','02e_z_9_6','02f_z_9_6','02g_z_9_6','02h_z_9_6','02i_z_9_6','02j_z_9_6']
    chains2=['02a_z_5x','02b_z_5x','02c_z_5x','02d_z_5x','02e_z_5x','02f_z_5x','02g_z_5x','02h_z_5x','02i_z_5x_2','02j_z_5x']
    chains2=['02a_z_6x','02b_z_6x','02c_z_6x','02d_z_6x','02e_z_6x','02f_z_6x','02g_z_6x','02h_z_6x','02i_z_6x','02j_z_6x' ]
    chains3='02z_31_2'
    print expt.parameters
    
    z_new=False
    Replot=False
    z_nov = [0.31, 0.5, 0.69, 0.9, 1.16, 1.44, 1.81, 2.18, 2.81, 3.71, 4.83]
    if z_new:
        z =[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,1.9,2.2,2.5, 2.8,3.2,3.6,4.0]
        chains=['11a','11b','11c_3','11d_3','11e_3','11f_3','11g_3','11h','11i','11j', '11k', '11l' ,'11m', '11n']#chains_191111
    else:
        z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]#oldest ;)
        
    sfrd, sfrd_l,sfrd_h,z_med = numpy.zeros(len(z)+1),numpy.zeros(len(z)+1) ,\
         numpy.zeros(len(z)+1) ,numpy.zeros(len(z)+1) 
    sfrd_nov, sfrd_nov_l,sfrd_nov_h, sfrd_2, sfrd_3, sfrd_2_l, sfrd_3_l, sfrd_2_h, sfrd_3_h  = numpy.zeros(len(z)+1),\
         numpy.zeros(len(z)+1), numpy.zeros(len(z)+1) ,numpy.zeros(len(z)+1),numpy.zeros(len(z)+1),\
         numpy.zeros(len(z)+1), numpy.zeros(len(z)+1) ,numpy.zeros(len(z)+1),numpy.zeros(len(z)+1) 
    
    
    ana3=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_2002%s'%chains3,outstem))
    drawmap3=ana3.get_best_fit()['parameters']
    params3 = drawmap3
    print 'chains_2002%s'%chains3
    print params3
    parameters2= ['noise', 'A_SF', 'A_agn', 'LMIN', 'LMAX']
    parameters3= ['noise', 'A_agn', 'A_SF', 'B_agn', 'B_SF', 'LMIN', 'LMAX']
    parameters3= ['noise', 'A_agn', 'A_SF', 'B_agn', 'B_SF', 'LMIN_1', 'LMIN_2', 'LMIN_3', 'LMIN_4', 'LMIN_5', 'LMIN_6', 'LMIN_7', 'LMIN_8', 'LMIN_9', 'LMIN', 'LMAX']
    for num in range(1,11):
    #num = 1.
     print num
     z_min,z_max, z_m = get_z(num,z_new)
     dl = get_dl(z_m)
     Vmax = get_Vmax(z_min,z_max)
     dsdl = get_dsdl(z_m,dl)

     z_m = (z_min + z_max)/2                       
     dl = get_dl(z_m)
    
     print z_min,z_m,z_max
     print expt.kind
     area = SURVEY_AREA*sqDeg2sr 
     area1 = SURVEY_AREA*sqDeg2sr
     print SURVEY_AREA
     ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_20%s%s'%(pre_chain,chains[num-1]),outstem))
     drawmap=ana.get_best_fit()['parameters']
     drawmid=ana.get_stats()['marginals']
     params = drawmap
     
     print params
     
     ana2=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_2002%s'%(chains2[num-1]),outstem))
     drawmap2=ana2.get_best_fit()['parameters']
     params2 = drawmap2
     
     print 'parameter 1 ', params
     print 'parameters 2 ', params2
     print 'parameters 3 ', params3
     
     
     
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
     
     if model in ['LFsch','LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFlognorm_dpl']:
        Lnorm=params[parameters.index('LNORM')]
        Lstar=params[parameters.index('LSTAR')]
        Lslope=params[parameters.index('LSLOPE')]
        #Lzevol=params[parameters.index('LZEVOL')]
     if model in ['LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFlognorm_dpl']:
        Lslope2=params[parameters.index('LSLOPE2')]
     if model in ['LFlognorm','LFpl', 'LFdpl', 'LFdpl_pl','LFpl_dpl', 'LFlognorm_dpl','LFdpl_dpl']:
        Lnorm_2=params[parameters.index('LNORM_2')]
        Lstar_2=params[parameters.index('LSTAR_2')]
        Lslope_2=params[parameters.index('LSLOPE_2')]

     if model in ['LFdpl_dpl','LFdpl']:
        Lslope2_2=params[parameters.index('LSLOPE2_2')]
    
     if model in ['LFlognorm','LFlognorm_dpl']:
        Lsigma = params[parameters.index('LSIGMA')]

     Lmin=params[parameters.index('LMIN')]
     Lmax=params[parameters.index('LMAX')]
     #Lmax2=params[parameters.index('LMAX2')]
     Lmax2=Lmax

    
     SMIN  = get_sbins(numpy.power(10,Lmin),z_m,dl)*1e6
     SMAX  = get_sbins(numpy.power(10,Lmax),z_m,dl)*1e6
     sigma,fsigma = numpy.log10(get_Lbins([SURVEY_NOISE,SURVEY_NOISE*5],z_m,dl,'muJy')*(1.4/3)**(-.7))
     sigma2,fsigma2 = numpy.log10(get_Lbins([450,2400],z_m,dl,'muJy'))
     print z_m,dl, sigma,fsigma
     print SURVEY_NOISE,SURVEY_NOISE*5,SURVEY_AREA
     
     s=numpy.loadtxt('chains_20%s%s/recon_stats.txt'%(pre_chain,chains[num -1]))
     xrecon=s[:-1,0]; yrecon=s[:-1,1]
     yrecon_d=s[:-1,2]; yrecon_u=s[:-1,3]
     yrecon_rms = s[:-1,4]
     yrecon_avr = s[:-1,5]
     yrecon_rms_down = yrecon_rms
     
     yrec_u = yrecon_u - yrecon
     yrec_d = yrecon - yrecon_d
     
     try:
        zrecon, srecon, srecon_l, srecon_h =numpy.loadtxt('sfrd_chains_2001%s_Lmin_%s'%(chains[0],'20'),unpack=True)
     except:
        print 'Please note that the 95% region has to be run... python reconstruct_lf_plot.py'
        print 'Ensure that chains in reconstruct_lf_plot.py cooresponds to chains in this file!'
        zrecon, srecon, srecon_l, srecon_h=0,0,0,0
     #print yrecon_d
     #print yrecon
     #print yrecon_u
     #sys.exit()
     
     #y_err_u = sqrt(yrec)     
     
     xreco = get_Lbins(xrecon,z_m,dl,'muJy')#*(1.4/3)**(-.7)
     xreco2 = get_Lbins(xrecon,z_m,dl,'muJy')
     yreco = yrecon
     lin_yrecon = numpy.log10(yreco)
     lin_xrecon = numpy.log10(xreco)
     lin_xrecon2 = numpy.log10(xreco2)
     q_l,q,q_h = get_q(z_m)
     
     #print z_nov[num -1]
     try:
        yrecon_nov = lumfuncUtils.lognormpl(get_Lbins(xrecon,z_m,dl,'muJy')/(1 + z_nov[num-1])**(3.16-z_nov[num-1]*0.32), numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
        
     except:
        print ''
     
     str_q_l = get_sfr_q(1.,q_l,xreco2)
     str_q_h = get_sfr_q(1.,q_h,xreco2)
     str_q   = get_sfr_q(1.,q,xreco2)
     
     
     #sys.exit()
     
     lmin = 20.#. 21.5
     Lmax=26.
     sfrd_z_l=get_sfrd_z(yrecon,lin_xrecon2,str_q_l,params,parameters,model,q_l ,Lmax2,lmin,z_m)[0] #lower limit
     sfrd_z_h=get_sfrd_z(yrecon,lin_xrecon2,str_q_h,params,parameters,model,q_h ,Lmax2,lmin,z_m)[0] #upper limit
     sfrd_z  =get_sfrd_z(yrecon,lin_xrecon2,str_q  ,params,parameters,model,q ,Lmax2,lmin,z_m)[0]  #26, 21.5
     
     sfrd_2_z =get_sfrd_z(yrecon,lin_xrecon2,str_q  ,params2,parameters2,'LFevol_logn_mat',q ,Lmax2,lmin,z_m)[0]
     sfrd_3_z =get_sfrd_z(yrecon,lin_xrecon2,str_q  ,params3,parameters3,'LFevol_logn',q ,Lmax2,lmin,z_m)[0]
     
     sfrd_2_z_l =get_sfrd_z(yrecon,lin_xrecon2,str_q_l  ,params2,parameters2,'LFevol_logn_mat',q_l ,Lmax2,lmin,z_m)[0]
     sfrd_3_z_l =get_sfrd_z(yrecon,lin_xrecon2,str_q_l  ,params3,parameters3,'LFevol_logn',q_l ,Lmax2,lmin,z_m)[0]
     sfrd_2_z_h =get_sfrd_z(yrecon,lin_xrecon2,str_q_h  ,params2,parameters2,'LFevol_logn_mat',q_h ,Lmax2,lmin,z_m)[0]
     sfrd_3_z_h =get_sfrd_z(yrecon,lin_xrecon2,str_q_h  ,params3,parameters3,'LFevol_logn',q_h ,Lmax2,lmin,z_m)[0]
     
     #sfrd_z_l=get_sfrd_z(yrecon,lin_xrecon2,str_q_l,params,parameters,model,q_l ,Lmax2,Lmin)[0] #lower limit
     #sfrd_z_h=get_sfrd_z(yrecon,lin_xrecon2,str_q_h,params,parameters,model,q_h ,Lmax2,Lmin)[0] #upper limit
     #sfrd_z  =get_sfrd_z(yrecon,lin_xrecon2,str_q  ,params,parameters,model,q ,Lmax2,Lmin)[0]  #26, 21.5
     #sys.exit()
     
     sfrd_nov_z =get_sfrd_z(yrecon_nov,lin_xrecon,str_q ,params,parameters,'novak',get_q(z_nov[num-1])[1] ,26,lmin, z=z_nov[num-1])[0] 
     print 'This is nov z ', z_nov[num-1]
     #sys.exit()
     #sfrd_z  =get_sfrd_z(yrecon,lin_xrecon,str_q  ,params,parameters,model,q ,Lmax2,Lmin)[0]
     
     
     #sfrd_z,sfrd_zerr= get_sfrd_z(yrecon,lin_xrecon,str_q,numpy.inf,0.)
     
     
     sfrd[num-1]   = sfrd_z
     sfrd_l[num-1] = sfrd_z_l
     sfrd_h[num-1] = sfrd_z_h
     z_med[num-1]    = z_m
     sfrd_nov[num-1] =sfrd_nov_z
     sfrd_2[num-1] = sfrd_2_z
     sfrd_3[num-1] = sfrd_3_z
     sfrd_2_l[num-1] = sfrd_2_z_l
     sfrd_3_l[num-1] = sfrd_3_z_l
     sfrd_2_h[num-1] = sfrd_2_z_h
     sfrd_3_h[num-1] = sfrd_3_z_h
     
     #plt.plot(lin_xrecon,numpy.log10(str_q),'*',markersize=8)

    
    #plt.rc('lines', linewidth=2)
    #plt.xlabel(r'$\rm{log_{10}[L_{1.4}/(W}$ $\rm{Hz^{-1})]}$',fontsize = 30)
    #plt.ylabel(r'$\rm{log_{10}[SFR_{IR}(M_\odot yr^{-1})]}$',fontsize = 30)
    #plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    #plt.tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    
    #plt.show()
    #fig = plt.figure()
    #print z_med
    #print sfrd_l
    #print sfrd
    #print sfrd_h
    #print len(z),len(sfrd)
    #Salpeter to Chabrier * 0.63
    
    z_b=numpy.arange(0,60)/10.
    Ilbert_h = sfrd_Behroozi(z_b,z0=0.95-0.410,A=-1.,B=0.194-0.082,C=0.111+0.04)
    Ilbert_l = sfrd_Behroozi(z_b,z0=0.95+0.343,A=-1.,B=0.194+0.128,C=0.111-0.029)
    G_x , Gruppioni       = numpy.loadtxt('Gruppioni_2013',unpack=True,delimiter=',')
    G_x1, Gruppioni_yup   = numpy.loadtxt('Gruppioni_2013_y_up',unpack=True,delimiter=',')
    G_x2, Gruppioni_ydown = numpy.loadtxt('Gruppioni_2013_y_down',unpack=True,delimiter=',')
    Burg_z, Burg_FIR, Burg_FIR_er, Burg_tot, Burg_tot_er = numpy.loadtxt('Burgarella_2013',unpack=True)
    
    
    plt.errorbar(z_med,numpy.log10(sfrd_2),fmt='--k',linewidth=5, label='Indvidual Novak LF Model fit')
    plt.errorbar(z_med,numpy.log10(sfrd_3),fmt='.-g',linewidth=6, label='Evolution Novak LF Model fit')
    try:
        plt.errorbar(numpy.array(nov_z),numpy.array(nov_sfrd),xerr=[nov_z_l,nov_z_u],yerr=[nov_sfrd_l,nov_sfrd_u], fmt='sr',fillstyle='none',ecolor='k',markeredgewidth=1.5, markersize=10, label='Novak et al. 2017 Pure L evoluion')
    except:
        print ''
    plt.errorbar(z_med,numpy.log10(sfrd),fmt='h',fillstyle='none',markeredgewidth=3,markersize=18, label='Our RLF model')
    #plt.errorbar(z_nov[:-1],numpy.log10(sfrd_nov[:-2]),fmt='sg',markersize=20, label='Novak et al. 2017 Pure L evoluion')
    #plt.errorbar(z_b,numpy.log10(sfrd_Behroozi(z_b)),fmt='s-',linewidth=2,fillstyle='none', label='Behroozi et al. 2013')
    #plt.errorbar(z_b,numpy.log10(sfrd_Behroozi(z_b,z0=0.95,A=-1.,B=0.194,C=0.111)),fmt='>:',linewidth=2,fillstyle='none', label='IIlbert et al. 2013 (UltraVista DR1)')
    #plt.errorbar(z_b,numpy.log10(sfrd_Madau(z_b)),fmt='o-',linewidth=2,fillstyle='none', label='Madau&Dickinson 2014')
    #plt.errorbar(z_med,numpy.log10(sfrd_l)+0.4,fmt='--',markersize=8,label='low')
    #plt.errorbar(z_med,numpy.log10(sfrd_h)+0.4,fmt='--',markersize=8,label='high')
   
    #plt.errorbar(numpy.array(nov_z),numpy.array(Lower),xerr=[nov_z_l,nov_z_u],yerr=[lw_up, lw_lw], fmt='*g', label='lower limit SFRD Novak et al 2017')
    #plt.yscale('log')
    #plt.fill_between(z_med,numpy.log10((sfrd_l)),numpy.log10(sfrd_h), color='b',alpha=0.3)
    plt.fill_between(z_med,numpy.log10((sfrd_2_l)),numpy.log10(sfrd_2_h), color='c',alpha=0.7)
    plt.fill_between(z_med,numpy.log10((sfrd_3_l)),numpy.log10(sfrd_3_h), color='g',alpha=0.5)
    plt.fill_between(zrecon,numpy.log10((srecon_l)),numpy.log10(srecon_h), color='b',alpha=0.3)
    #plt.fill_between(z_b,numpy.log10((Ilbert_l)),numpy.log10(Ilbert_h), color='orange',alpha=0.2)
    plt.fill_between(G_x,numpy.log10(Gruppioni_ydown),numpy.log10(Gruppioni_yup), color='k',alpha=0.2,label='Gruppioni-2013 IR')
    plt.fill_between(nov_z,nov_phi_L, nov_phi_L2, color='r',alpha=0.4,label='Novak et al. 2017 $\Phi$ and L evolution')
    #plt.fill_between(Burg_z,numpy.log10(0.63*0.01*(Burg_FIR - Burg_FIR_er)),numpy.log10(0.63*0.01*(Burg_FIR + Burg_FIR_er)), color='r',alpha=0.2,label='Burgarella-2013 FIR')
    #plt.fill_between(Burg_z,numpy.log10(0.63*0.01*(Burg_tot - Burg_tot_er)),numpy.log10(0.63*0.01*(Burg_tot + Burg_tot_er)), color='g',alpha=0.2,label='Burgarella-2013 Total')
    plt.ylabel(r'$\rm{\log10(SFRD[M_\odot yr^{-1} Mpc^{-3}])}$',fontsize=30)
    plt.xlabel(r'$\rm{z} $',fontsize = 30)    
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
    plt.legend(frameon=False, prop={'size':20}).draggable()
    plt.show()
    
	
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
