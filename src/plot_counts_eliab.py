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
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats,calculate_confidence2,calculate_confidence,peak_confidence
from countUtils import calculateDnByDs,medianArray
from lumfuncUtils import get_z, get_dl, get_sbins
from plat import open_anytype

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

#-------------------------------------------------------------------------------
def dnds(flux,area,bins=numpy.logspace(-6,1,30),error=False):
    '''
    Computes the eucledian weighted dnds
    parameters; flux area and bins
    
    flux: the flux in uJy
    area: survey area in str
    bins: flux bins in uJy
    
    returns
    x: Medianbins in uJy
    y: dnds*s^2.5
    
    '''
    print bins
    h = numpy.histogram(flux,bins)
    print h
    x = medianArray(h[1]*1e6)
    y = calculateDnByDs(h[1],h[0],eucl=False,idl_style=False) /area 
    if error:
     err =calculateDnByDs(h[1],h[0],errors=True)/area
     return x,y,err
    return x,y
def multiply(a,b):
    N=len(a)
    c=numpy.zeros(N)
    for i in range(N):
        c[i] = a[i]*b[i]
    return c
#-------------------------------------------------------------------------------
@profile
def main():

    """
    """

    print 'Settings file is %s' % setf

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    bins2 = [ -6.53657634e+02,  -2.60225791e+02,  -1.03597753e+02,\
        -4.12430085e+01,  -1.64191374e+01,  -6.53657634e+00,\
        -2.60225791e+00,  -1.03597753e+00,  -4.12430085e-01,\
        -1.64191374e-01,  -6.53657634e-02,  -2.60225791e-02,\
        -1.03597753e-02,  -4.12430085e-03,  -1.64191374e-03,\
         1.64191374e-03,   4.12430085e-03,   1.03597753e-02,\
         2.60225791e-02,   6.53657634e-02,   1.64191374e-01,\
         4.12430085e-01,   1.03597753e+00,   2.60225791e+00,\
         6.53657634e+00,   1.64191374e+01,   4.12430085e+01,\
         1.03597753e+02,   2.60225791e+02,   6.53657634e+02,\
         1.64191374e+03,   4.12430085e+03,   1.03597753e+04,\
         2.60225791e+04,   6.53657634e+04,   1.64191374e+05,\
         4.12430085e+05,   1.03597753e+06,   2.60225791e+06,\
         6.53657634e+06]
    bins1 = numpy.arange(15.,29.2,0.4)
    print 'run expt'
    nbins = len(bins1)

    z_m = redshifts[0]
    dl = get_dl(z_m) 
    #z_min,z_max, z_m = get_z(num)    
    sbin1 = get_sbins(10**bins1,z_m,dl)*1e6
    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise,doRedshiftSlices=True,mybins=sbin1)
    expt2=countModel(modelFamily,nlaws,setf,[dataset],floatNoise,doRedshiftSlices=True)#,mybins=sbin1)
    #print expt.data
    print expt2.bins
    print expt2.data
    #sys.exit()

    # Get MAP parameters
    print '\n'*2
    print len(expt.parameters)
    print expt.parameters
    ncols = len(expt.parameters)

    f='%spost_equal_weights.dat' % outstem
    f=os.path.join(outdir,f)
     # Fetch best-fit parameters and calculate best-fit line
    
    # Load equally-weighted posterior samples
    x=numpy.genfromtxt(f)
    nsamp=x.shape[0]
    ncols=x.shape[1] 
    z=numpy.zeros((nsamp,ncols-1+expt.nbins))
    z[:,:-(expt.nbins-1)]=x
    # Shift posterior values to end
    z[:,-1]=z[:,ncols-1] # Copy...f='%spost_equal_weights.dat' % outstem
    z[:,ncols-1]=0.0 

     # Fetch best-fit
    summf=os.path.join(outdir,'1-summary.txt')
    summary=numpy.genfromtxt(summf)[-1,:]
    drawmap=summary[-(ncols+1):-2] 
    
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawml=ana.get_best_fit()['parameters']
    drawmap=drawml
    print drawmap
    Lmin=10**drawmap[1]
    Lmax=10**drawmap[2]
    [SMIN,SMAX]=get_sbins([Lmin,Lmax],z_m,dl)

    ymap=expt.evaluate2(drawmap)[0]
    #yram=expt.evaluate2([150.,19.,25.,-2.6,22.6,1.28,0.02])[0]
    
    #s25 = numpy.power(expt.binsMedian/1.0e6,2.5)
    #ymap= multiply(ymap,s25)
    print ymap[10]
    print 'good job'
    #print drawmap
    
    print nsamp,ncols, len(expt.bins)
    for isamp in xrange(nsamp):
        #z[isamp,ncols-1:]=expt.evaluate(expt.convertPosterior(z[isamp,:],power))
        z[isamp,ncols-1:]=expt.evaluate2(z[isamp,:])[0]#,expt.binsMedian)
        #z[isamp,ncols-1:]=multiply(z[isamp,ncols-1:], s25)
        #print z[isamp,:]
        #sys.exit()

    # Blanking, 0.0 -> NaN
    z[numpy.where(z==0.0)]='NaN'
    recons=z[:,ncols-1:]
    #print recons[0]
    print 'after'
    print z[0]
    print 'recon'
    print recons[0]
    #sys.exit()

    s=numpy.zeros((len(expt.binsMedian),8))
    s[:,0]=expt.binsMedian

    print '# ibin flux fit low high dlower dupper skew kurtosis'
    #print expt.binsMedian
    #sys.exit()
    for ibin in xrange(len(expt.binsMedian)):
        x = recons[:,ibin]
        # Remove NaNs from stats vectors
        # http://stackoverflow.com/questions/11620914/removing-nan-values-from-an-array
        x = x[~numpy.isnan(x)]
        #ss=stats.bayes_mvs(x,alpha=0.68)[0]
        #x*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print x
        #sys.exit()

        try:
            ss=numpy.zeros(3)
            #print 'worked',ibin,expt.binsMedian[ibin]
            ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence(x,alpha=0.95,ret_all=True)
            #ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence2(x,alpha=0.95,ret_all=True,\
            #                                                   value_central=ymap[ibin],\
            #                                                   truncate_edges=True)
        except:
            ss=numpy.nan*numpy.ones(3)
            print 'not',expt.binsMedian[ibin]
            #print "didn't work for ", x
            continue
        #sys.exit()
        tt=peak_confidence(x,bins=10)
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
#        s[ibin,1]=ss[0]  # median
#        s[ibin,1]=tt     # peak
        s[ibin,1]=ymap[ibin] # MAP
        s[ibin,2]=ss[1]  # lower
        s[ibin,3]=ss[2]  # upper
        s[ibin,4]= 2*numpy.std(x)# (s[ibin,3]- s[ibin,2])/(2*numpy.std(x))# skewness
        s[ibin,5]= ss[0]
        
    
    if dataset=='sdss':
        print datafiles, len(datafiles[0])
        dataf = datafiles[0]
        if len(dataf)==20:
            num = dataf[-5]
            print num
        else:
            num = dataf[-6:-4]
        d=numpy.genfromtxt('%s/%s'%(outdir,datafiles[0][4:]))
        bins=d[:,2]
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
        bins=d[:,2]
        print datafiles[0][8:]

    elif 'sim' in dataset:
        d=numpy.genfromtxt(os.path.join(outdir,'sim.txt'))
        #dnds_d = d[:,5]
    num = int(num)
    #num = 1.
    print num
    z_min,z_max, z_m = get_z(num)    
    
    SURVEY_AREA=1.38
    #flux2 = open_anytype('../First_dr12d_1_8arc.txt',(2,14),[z_min,z_max],S=True)*1e-3
    #flux2 =open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,21),[z_min,z_max],F='uJy',S=True)
    #flux3 =open_anytype('cos_data/cos_s1.txt',(2,3),[z_min,z_max],F='mJy',S=True)
    flux3 = open_anytype('cos_data/cos_s%s.txt'%num,(2,3), F='uJy',S=True,band='S')
    
    print min(bins)*1e-6,max(bins)*1e-6    
    nbin = 20
    #sys.exit()

    #bins =numpy.array(bins)
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_d=s[:-1,2]; yrecon_u=s[:-1,3]
    yrecon_rms = s[:-1,4] 
    yrecon_avr = s[:-1,5] 
    yrecon_rms_down = yrecon_rms
    
    
    print type(flux3)
    print flux3
    
    #sys.exit()
    #x4,y4,err4 = dnds(flux2,SURVEY_AREA*sqDeg2sr,bins=expt2.binsMedian[expt2.binsMedian>0]*1e-6,error=True)
    #x5,y5,err5 = dnds(flux3,SURVEY_AREA*sqDeg2sr,bins=expt2.binsMedian[expt2.binsMedian>0]*1e-6,error=True)
    #x4,y4,err4 = dnds(flux2*1e-6,SURVEY_AREA*sqDeg2sr,bins=expt2.bins[expt2.bins>0]*1e-6,error=True)
    #x4,y4,err4 = dnds(flux2*1e-6,SURVEY_AREA*sqDeg2sr,bins=sbin1*1e-6,error=True)
    x5,y5,err5 = dnds(flux3*1e-6,SURVEY_AREA*sqDeg2sr,bins=sbin1*1e-6,error=True)
    #x6,y6,err6 = dnds(flux4*1e-6,SURVEY_AREA*sqDeg2sr,bins=sbin1*1e-6,error=True)
    #x5,y5,err5 = dnds(flux3*1e-6,SURVEY_AREA*sqDeg2sr,bins=expt2.bins[expt2.bins>0]*1e-6,error=True)
    #fig = plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    
    dn_by_ds,dn_by_ds_eucl,dn_by_ds_errs,dnbdsb,j,ds =expt2.dn_by_ds(return_all=True,data=expt2.data) 
    # Plot the data
    #print expt.binsMedian
    #print dn_by_ds_eucl
    #print expt.data
    #plt.plot(x0,y0,'rx',label='Extracted fluxes')    
    print len(expt.bins),len(dn_by_ds_eucl)
    print len(ymap),len(xrecon),len(expt.bins), len(expt.binsMedian)
    #plt.errorbar(x4,y4,yerr=err4,fmt='sb',label='Detected',markersize=7)
    #plt.plot(expt.binsMedian,numpy.array(yram)/SURVEY_AREA/sqDeg2sr,'or',label='by-eye fit')
    #plt.errorbar(x6,y6,yerr=err6,fmt='sr',label='Detected > 1muJy',markersize=7)
    plt.errorbar(x5,y5,yerr=err5,fmt='ok',label='Extracted',markersize=6)
    plt.errorbar(ds*1e6,dn_by_ds/SURVEY_AREA/sqDeg2sr,yerr=dn_by_ds_errs,fmt='*k',label='extracted',markersize=10)
    plt.plot(xrecon,yrecon/SURVEY_AREA/sqDeg2sr,'--k', label='MAP',linewidth=2)
    plt.plot(xrecon,yrecon_avr/SURVEY_AREA/sqDeg2sr,'--c', label='Average',linewidth=2)

    plt.fill_between(xrecon[xrecon>SMIN],\
                     (yrecon_d[xrecon>SMIN])/SURVEY_AREA/sqDeg2sr\
                     ,yrecon_u[xrecon>SMIN]/SURVEY_AREA/sqDeg2sr,\
                     color='b',alpha=0.3)
                     
    #plt.fill_between(xrecon[xrecon>SMIN],\
    #                 (yrecon_avr[xrecon>SMIN]-yrecon_rms[xrecon>SMIN])/SURVEY_AREA/sqDeg2sr\
    #                 ,(yrecon_avr[xrecon>SMIN]+\
    #                 yrecon_rms[xrecon>SMIN])/SURVEY_AREA/sqDeg2sr,\
    #                 color='r',alpha=0.3)                 
    
    plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
   # 		plt.text(8e3,3e-6,models[i],fontsize=24)
    plt.text(2e3,1e-3,'%.2f < z < %.2f'%(z_min,z_max),fontsize=16,alpha=0.7 )
    f=30
    #legend=plt.legend(loc='lower right',prop={'size':15},frameon=False,numpoints=1)
    legend=plt.legend(loc='upper right',prop={'size':15},frameon=False,numpoints=1)
    #frame = legend.get_frame()
    #frame.set_facecolor('white')
    plt.xlim(1e-1,5e4)
    plt.ylim(7e4,5e13)
    plt.xlabel('$S(\mu \mathrm{Jy}^{-1})$',fontsize=f)
    plt.ylabel('$n(\mathrm{sr}^{-1}\mathrm{Jy}^{-1})$',fontsize=f)
    plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    plt.tick_params(axis='both',which = 'minor', labelsize=18, width=2)
    #plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()


    plotf='%s/counts.png' % (outdir)
    plt.show()
    #plt.savefig(plotf,bbox_inches='tight')
    print '-> Look in open %s' % plotf 

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
