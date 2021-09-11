"""
Support functions for bayestack, bayestackClasses

Luminosity functions and evolution thereof

No longer depends on dnds_lumfunc module

Jonathan Zwart
November 2015

"""

import os,sys
import importlib
import glob
import numpy
import mpmath
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.interpolate import interp1d
from scipy.special import erf
import matplotlib.pyplot as plt
from profile_support import profile
from utils import sqDeg2sr,sqrtTwo,find_nearest,medianArray,\
                           interpol,buildCDF,Jy2muJy,interpola
from cosmocalc import cosmocalc

if 'chains' in sys.argv[-1]:
    #potential_settings=glob.glob(os.path.join(sys.argv[-1][:sys.argv[-1].index('/')],'*settings*py'))
    potential_settings=glob.glob(os.path.join(sys.argv[-1],'*settings*py'))
    #print 'setting file ',potential_settings[0][:-3]
    assert(len(potential_settings)==1), '***More than one potential settings file!'
    #settingsf=potential_settings[0][:-3]
    settingsf='.'.join([sys.argv[-1],potential_settings[0].split('/')[-1].split('.')[-2]])
else:
    settingsf=sys.argv[-1].split('.')[-2]

#print '%s is using %s' % (__name__,settingsf)
try:
    set_module=importlib.import_module(settingsf)
    print settingsf
    globals().update(set_module.__dict__)
except:
    print '***Warning: Settings not loaded'

#-------------------------------------------------------------------------------
# LF utility functions

@profile
def erfss(S,Sbinlow,Sbinhigh,ssigma):
    """
    For a Gaussian noise model
    """
    if dataset=='sdss':
    	S = float(max(S/1.4,S-2.5e-4))
    else:
	S = float(S)
    return 0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

#-------------------------------------------------------------------------------

def get_Lbins(sbins,z,dl):
    """
    Units:
        sbins - Jy
        dl - Mpc
        Lbins - W Hz^-1
    """
    if isinstance(z,float):
       Lbins = [math.pi*4 *(s*1e-26)* (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1) for s in sbins]
    else:
        Lbins = [math.pi*4 *(s*1e-26)* (dli* 3.08e22)**2 * (1 + zi)**((LF_SPEC_INDEX)+1) for s,zi,dli in zip(sbins,z,dl)]
    return Lbins #W/Hz


#-------------------------------------------------------------------------------

def get_z(ind):
    """
Computes the redshift slice from the binfile

returns z_min, z_max and z_mean
    """
    z = [0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
    z_mid =[0.32, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
    z_min  = z[ind -1]
    z_max  = z[ind]
    z_mean = numpy.mean((z_min,z_max))
    z_mean = z_mid[ind -1]
    return z_min,z_max, z_mean
def get_num(z_i):	
   '''
   returns the index number corrisponding to z
   '''
   z_i=round(z_i,2)
   z = [0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
   z_mid =[0.32, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
   return z_mid.index(z_i)


#-------------------------------------------------------------------------------

def get_dl(z):
    """
Returns the comoving radial distance in Mpc
    """
    if isinstance(z,float):
       return cosmocalc(z,H0=Ho,WM = wm)['DCMR_Mpc']       
    else:
        
       return [cosmocalc(zi,H0=Ho,WM = wm)['DCMR_Mpc'] for zi in z]
    
#-------------------------------------------------------------------------------

def get_dL(z):
    """
Returns the comoving radial distance in Mpc
    """
    if isinstance(z,float):
       return cosmocalc(z,H0=Ho,WM = wm)['DL_Mpc']       
    else:
        
       return [cosmocalc(zi,H0=Ho,WM = wm)['DL_Mpc'] for zi in z]
    
#-------------------------------------------------------------------------------


def get_dsdl(z,dl):
    """
    dl is luminosity distance in Mpc
    dsdl is dS/dL in metres^-2
    """
    dsdl = 1./(math.pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    print 'dsdl'
    print dsdl
    return dsdl
#-------------------------------------------------------------------------------


    
def get_sbins(fbins,z,dl):
    """
    Units:
        fbins - W Hz^-1 [luminosity bins]
        dl - Mpc
        sbins - ***Jy
    """
    fbins = numpy.array(fbins)
    sbins = fbins/(pi*4 * (dl*3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX+1))) #sbins = fbins/(pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    return sbins*1e26 #1e26 #*10.9837 #Jy


#-------------------------------------------------------------------------------

def get_dlds(z,dl):
    """
    dl is comoving distance in Mpc
    dlds is dL/dS in metres^-2
    """
    dlds = (pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    #print dlds
    return dlds

#-------------------------------------------------------------------------------

def get_Vmax(zlo,zup):
    z  = zup
    V1 = cosmocalc(z,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
    V2 = cosmocalc(zlo,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
    area = SURVEY_AREA*sqDeg2sr
    vmax = area*(V1-V2)/(4*pi)
    return vmax

#-------------------------------------------------------------------------------

def LFtodnds(lbins, LF, dlds, V, dl):
    """
    Units:
        lbins   - W Hz^-1
        LF(phi) - Mpc^-3 mag^-1
        dlds    - Jy W^-1 Hz
        V(Vmax) - Mpc^-3
        dl      - Mpc
    Returns:
        dnds - 
    """

    L = numpy.array(lbins)
    phi_m = LF/(numpy.log(10)*L) #dLogL
    dndl = phi_m * V
    dnds = dndl*dlds
    return dnds*1e-26

#-------------------------------------------------------------------------------

@profile
def schechter(L,ln_Lstar,alpha,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha - no units
        * norm - phi_star Mpc^-3 mag^-1
    Outputs:
        phi
    """
    
    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star*over_Lstar
    phi = norm *(Lr**alpha) *numpy.exp(-Lr)
    return phi

#-------------------------------------------------------------------------------
def lognormpl(L,ln_Lstar,alpha,sigma,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha - no units
        * norm - phi_star Mpc^-3 mag^-1
        * sigma - no units
    Outputs:
        phi
    """
    
    #print 'parameters'
    
    
    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    sig_1= 1./(2*sigma**2)
    phi = norm *Lr**(1-alpha) * numpy.exp(-sig_1 *numpy.log10(1+ Lr)**2)
    return phi

@profile
def doublepowerlaw(L,ln_Lstar,alpha1,alpha2,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * alpha2 - no units
        * C - phi_0 Mpc^-3 mag^-1
    Outputs:
        phi - Mpc^-3 mag^-1
    """
    
    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    phi = norm* ((Lr)**alpha1 + (Lr)**alpha2)**-1.
    return phi

@profile    
def powerlaw(L,ln_Lstar,alpha1,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * C - phi_0 Mpc^-3 mag^-1
    Outputs:
        phi - Mpc^-3 mag^-1
    """

    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    phi = norm* ((Lr)**alpha1)**-1.
    return phi

    
@profile
def doublepower2law(L,Lstar,alpha1,alpha2,phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * alpha2 - no units
        * C - phi_0 Mpc^-3 mag^-1
    Outputs:
        phi - Mpc^-3 mag^-1
    """
    
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    phi = norm* ((Lr)**alpha1 + (Lr)**alpha2)**-1.
    return phi
    
def ddoublepowerlaw(L,ln_phi_star,ln_phi_star2,alpha,ln_L1,beta,ln_L2,gamma,ln_L3,delta):
    '''
    Double double power law
     Inputs:
        Lbins  - Log10(W Hz^-1)
        *L1    - Lstar1,first break [W Hz^-1]
        *L2    - second break 
        *L3    - Lstar2, third break [W Hz^-1]
        *alpha - Slope  
        *beta,gamma,delta - second,third and fouth slope [no units]
        * phi_star - normalization [phi_0 Mpc^-3 mag^-1]
    Outputs:
        phi - Mpc^-3 mag^-1

    '''
    phi_star = numpy.power(10,ln_phi_star)
    phi_star2 = numpy.power(10,ln_phi_star2)
    Lstar1    = numpy.power(10,ln_L1)
    Lstar2    = numpy.power(10,ln_L3)
    L2        = numpy.power(10,ln_L2)  
    
    if L <=L2:
        Lr    = L/Lstar1
        phi   = phi_star* ((Lr)**alpha + (Lr)**beta)**-1.
    elif L2 < L:
        Lr    = L/Lstar2
        phi   = phi_star2* ((Lr)**gamma + (Lr)**delta)**-1.
        
    return phi
#-------------------------------------------------------------------------------

@profile
def dNdS_LF(S, z, dlds, Vmax, dl, params=None, paramsList=None, family=None):
    """
    Source count model
    S is in Jy
    Cosmology is set globally
    """
    
    L=get_Lbins([S],z,dl)[0]*(1.4/3.)**(-0.7) #doing everything at 1.4GHz
    phi = LF(S, z, dlds,Vmax,dl,params,paramsList,family,L=L)
    return LFtodnds(L*(1.4/3.)**(0.7),phi,dlds,Vmax,dl)
   
def LF(S, z, dlds, Vmax, dl, params=None, paramsList=None, family=None,L=None):
    '''
    Calculates the LF given the model, parameters and redshift.
    '''
    #setting up the params
    Lmin=params[paramsList.index('LMIN')]
    Lmax=params[paramsList.index('LMAX')]
    if family in ['LFsch','LFdpl_dpl','LFlognorm_dpl','LFlognorm_lognorm_dpl']:
        Lnorm=params[paramsList.index('LNORM')]
        Lstar=params[paramsList.index('LSTAR')]
        Lslope=params[paramsList.index('LSLOPE')]
        #Lzevol=params[paramsList.index('LZEVOL')]
    if family in ['LFdpl_dpl','LFlognorm_dpl','LFlognorm_lognorm_dpl']:
        Lslope2=params[paramsList.index('LSLOPE2')]
    if family in ['LFlognorm','LFpl', 'LFlognorm_dpl','LFdpl_dpl','LFlognorm_lognorm_dpl']:
        Lnorm_2=params[paramsList.index('LNORM_2')]
        Lstar_2=params[paramsList.index('LSTAR_2')]
        Lslope_2=params[paramsList.index('LSLOPE_2')]
    if family in ['LFlognorm_lognorm_dpl']:
        Lnorm_3=params[paramsList.index('LNORM_3')]
        Lstar_3=params[paramsList.index('LSTAR_3')]
        Lslope_3=params[paramsList.index('LSLOPE_3')]

    if family in ['LFlognorm_dpl','LFlognorm_lognorm_dpl','LFdpl_dpl']:
        Lmin2=params[paramsList.index('LMIN2')]
        Lmax2=params[paramsList.index('LMAX2')]
        Lmin2,Lmax2 =numpy.power(10,[Lmin2,Lmax2])

    if family in ['LFdpl', 'LFdpl_dpl']:
        Lslope2_2=params[paramsList.index('LSLOPE2_2')]

    if family in ['LFevol_logn','LFevol_logn_L']:
	    alpha_agn=params[paramsList.index('A_agn')]
	    alpha_SF=params[paramsList.index('A_SF')]
	    beta_agn=params[paramsList.index('B_agn')]
	    beta_SF=params[paramsList.index('B_SF')]

    if family in ['LFevol_dpl_s','LFevol_logn_s']:
        alpha_SF=params[paramsList.index('A_SF')]
        beta_SF=params[paramsList.index('B_SF')]


    if family in ['LFevol_logn_el']:
        alpha_SF=params[paramsList.index('A_SF')]
        alpha_agn=params[paramsList.index('A_agn')]

    if family in ['LFevol_logn_lmin']:
	    Lmin2=params[paramsList.index('LMIN2')]
    

    if family in ['LFevol_logn_L']:
	    num=get_num(z)
        if num>0:
	        Lmin=params[paramsList.index('LMIN_%d'%num)]
	    
    if family in ['LFlognorm','LFlognorm_dpl','LFpl_lognorm','LFevol_logn_all','LFevol_logn_all_L']:
        Lsigma = params[paramsList.index('LSIGMA')]
    #print params
    if L==None:
    	L=get_Lbins([S],z,dl)[0]*(1.4/3.)**(-.7)        
    
    Lmin,Lmax =numpy.power(10,[Lmin,Lmax])
    
    if Lmin < L < Lmax:
            #print  Lmin,Lmax, L
            if family=='LFsch':
                phi=schechter(L,Lstar,Lslope,Lnorm)
            elif family=='LFpl':
                phi=powerlaw(L,Lstar_2,Lslope_2,Lnorm_2)            
            elif family=='LFdpl':
                if  Lslope_2 < Lslope2_2:
                    return -1.0e99
                phi=doublepowerlaw(L,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
            elif family=='LFdpl_pl':
                if  Lslope < Lslope2:
                    return -1.0e99

                phi_1,phi_2 =0,0
                if Lmin2<L<Lmax:
                    phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                if Lmin<L<Lmax2:
                    phi_2=powerlaw(L,Lstar_2,Lslope_2,Lnorm_2)
                phi = phi_1+phi_2


            elif family=='LFdpl_dpl':
                #Model A
		        phi_1,phi_2 =0,0
		        if Lmin2<L<Lmax:
		            phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
		        if Lmin<L<Lmax2:
			    phi_2=doublepowerlaw(L,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
		        phi = phi_1+phi_2


            elif family in ['LFevol_logn', 'LFevol_logn_L'] :
                        #Model C PLE
                        L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
                        phi = phi_1+ phi_2

                        phi = phi_1+phi_2
	   
	    elif family=='LFevol_logn_el':
	            #Model C indivi
                L_1 = L/(1+z)**(alpha_agn)
                L_2 = L/(1+z)**(alpha_SF)
                phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
			    phi_2=lognormpl(L_2,21.36,1.15,0.5, -2.36) #z_6f
			    phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                phi = phi_1+phi_2
	

            elif family=='LFlognorm_lognorm_dpl':
                #model D
	            phi_1,phi_2,phi_3 =0,0,0
       	       	if Lmin2<L<Lmax:
                    phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                if Lmin<L<Lmax2:
                    phi_2=lognormpl(L,Lstar_2,Lslope_2,0.6,Lnorm_2)
                    phi_3=lognormpl(L,Lstar_3,Lslope_3,0.6,Lnorm_3)
                phi = phi_1+phi_2+phi_3


            elif family=='LFlognorm':
                phi=lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
            elif family=='LFlognorm_dpl':
                #Model B
	            phi_1,phi_2 =0,0
       	       	if Lmin2<L<Lmax:
                    phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                if Lmin<L<Lmax2:
                    phi_2=lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
                phi = phi_1+phi_2

            return phi
    else:
    	#print 'passed', Lmin,L,Lmax
        return 0.

#-------------------------------------------------------------------------------

@profile
def IL(dNdS_LF,redshift,dlds,Vmax,dl,params,paramsList,Sbinlow,Sbinhigh,\
       inta=None,area=None,family=None,NOISE=None):
    """
    The integrand is the product of dN/dS_LF (count model) and G (the Gaussian)
    Schechter: ['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL'] + ['noise']
    Double PL: ['LNORM','LSTAR','LSLOPE','LSLOPE2','LMIN','LMAX','LZEVOL'] + ['noise']
    """
    
    Lmin=10**params[paramsList.index('LMIN')]
    Lmax=10**params[paramsList.index('LMAX')]
    sigma=params[paramsList.index('noise')]
    if family in ['LFevol_dpl_L','LFevol_logn_L']:
        num=get_num(redshift)
        if num>0:
          Lmin=10**params[paramsList.index('LMIN_%d'%num)]
    [Smin,Smax]=get_sbins([Lmin,Lmax],redshift,dl)

    #print redshift,dl,Sbinlow,Sbinhigh,Lmin,Lmax,Smin,Smax
    #intt = integrate.quad(lambda S:dNdS_LF(S,redshift,dlds,Vmax,dl,params=params,paramsList=paramsList,\
    #         family=family)*erfss(S,Sbinlow,Sbinhigh,sigma/1.0e6),\
    #                                    Smin,Smax,epsabs=1.49e-10, epsrel=1.49e-10,limit=600)[0]

    #I am intergrating in log space as it much faster. I've ran a may tests and varified that the two approaches give the same results.
    intt = integrate.quad(lambda x:dNdS_LF(exp(x),redshift,dlds,Vmax,dl,params=params,paramsList=paramsList,\
             family=family)*erfss(exp(x),Sbinlow,Sbinhigh,sigma/1.0e6)*exp(x),\
                                        log(Smin),log(Smax))[0]
    return intt
#-------------------------------------------------------------------------------

@profile
def calculateL3(params,paramsList,bins=None,\
                family=None,dump=None,verbose=False,dlds = None,redshift=None,Vmax=None,dl=None,noise=0):

    """
    For LF,
    function to calculate mock data for a given power law or interpolation object
    """

    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        IIi = abs(IL(dNdS_LF,redshift,dlds,Vmax,dl,params,paramsList,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,\
                    family=family,NOISE=noise))
        #II[ibin]=abs(IL(dNdS_LF,redshift,dlds,Vmax,dl,params,paramsList,bins[ibin],bins[ibin+1],\
        #            family=family))
        II[ibin] = IIi

    return II

#-------------------------------------------------------------------------------
