"""
Parameter file and inits for bayestack.py, lumfunc.py and simulate.py
"""

import os,sys,numpy
from utils import sqDeg2sr,beamFac

#-------------------------------------------------------------------------------
# Which program is loading this settings file?
exe=sys.argv[0]
context='cmdline'
if    'reconstruct' in exe: context='r'
elif  'simulate'    in exe: context='s'
elif  'lumfunc'     in exe: context='l'
elif  'bayestack'   in exe: context='y'
elif  'plot'        in exe: context='p'
elif  'binner'      in exe: context='b'
elif  'extractor'   in exe: context='e'
elif  'resample'    in exe: context='u'
elif  'inject'      in exe: context='i'
elif  'dnds_lumfunc'in exe: context='z'
elif  'picker'      in exe: context='w'
print 'Context is %s' % context
#-------------------------------------------------------------------------------

# New-style settings <-- bayestack.py

#dataset='video'
binStyle=1
nlaws=3 # Number of ppl laws (or poly coefficients) to fit
floatNoise=False
modelFamily='LFdpl'#'ppl' 'poly' # 'LFsch' # 'LFdpl'
outdir='chains_151109a' # based on 140123a
noiseZones=[0] # Select which noise zones to analyse
doPoln=False
doRayleigh=True # Set Rice signal to zero, i.e. Rayleigh distribution
doRedshiftSlices=True # For joint z slice analysis
whichRedshiftSlices=[1]
setBreaks=False
# LF variables
LF_SPEC_INDEX=-0.7
Ho=71.0
wm=0.27

simFamily= 'poly' # 'array' 'skads' 'ppl' 'poly' 'bins' 'test'
SMIN_SIM=0.5#0.01 # uJy
SMAX_SIM=85.0 # uJy
#simParams=[SMIN_SIM,SMAX_SIM]
#simParams=[1000.0,SMIN_SIM,SMAX_SIM,0.0]
#simParams=[SMIN_SIM,SMAX_SIM,1.0]
simParams=[SMIN_SIM,SMAX_SIM,5.0,1.0,0.5,0.1]
#simParamsList=['S0','S1']
#simParamsList=['C','S0','S1','a0']
#simParamsList=['S0','S1','p0']
simParamsList=['S0','S1','p0','p1','p2','p4']
simBins=numpy.linspace(-40.0,85.0,21)
SEED_SIM=1234
NSIM=72000
NOISE_SIM=16.2 # uJy; NB This is actually reset below...
dump='R.txt'
output='dummy.txt'
verbose=True
skadsFile='skads/1sqdeg_0p02uJy.txt'
simArrayFile='sims/150625a/sim_noiseless.txt'
simPolePosns=None
#SIM_DO_CAT_NOISE=True

#-------------------------------------------------------------------------------

# Master parameters
#MOTD=''
RESUME=False # Turn checkpointing on
#nb= 22#40#38#40#39#37#41#50#39#37 #13 #24 #27 #34 #37  # 38 or 41
dnds0=False # Leave as False otherwise no recon line...
binsHigh=False # Run with this True then set to False
#outdir='chains_150514a' # based on 140123a
run_num=outdir.split('_')[-1]
#if context=='l': outdir=os.path.join('/home/jtlz2/bartolomeu/output',outdir)
if context=='s' or context=='i': outdir='sims/%s' % outdir.split('_')[-1]
logfile='README.txt'
variablesfile='variables.txt'
comment='SDSS QSO RLF - tests'
SEED_SAMP=1234 # [-1 for clock]
#floatNoise=False

# What to fit - SPL, TPL, XPL or QPL
#nlaws=1

# Data set
#dataset='sims/151028a'
#dataset='cosmos'
#dataset='vvdf'
#dataset='video'
#run_num_run='151028a'
#dataset='first'
#dataset='mca'
#dataset='en1jvla'
dataset='sdss'

#-------------------------------------------------------------------------------
# Extraction settings
perGalMode=False
if True or context=='e':
    perGalMasterCat='/Users/jtlz2/video/20131007/VIDEO-cfhtls-D1_2013-10-07_fullcat_spec_fluxes_w_radio_jz_for_bayestack_140526.dat'
    # Select sample
    perGalGalaxyType='all'
#    perGalGalaxyType='ell'
#    perGalGalaxyType='mid'
#    perGalGalaxyType='sbn'

    perGalRadioPosn=False # Use radio or optical posn for extraction
    perGalWrapRoot='video/%s' % perGalGalaxyType
    #perGalCatForm=0
    #perGalCat='cats/smap_scat.txt'

    useInjectedMap=False
    #perGalUpsampleRate=None # if not None, do this on-the-fly; -ve -> read from file
    perGalUpsampleRate=-8 # if not None, do this on-the-fly; -ve -> read from file
    if perGalUpsampleRate is not None:
        injectedMap='sim_map_x%s.fits' % abs(perGalUpsampleRate)
    else:
        injectedMap='sim_map.fits'
    injectedMapDir=os.path.join('/Users/jtlz2/video/lumfunc-archive/sims-archive',dataset)
    positionsf='posns.txt'
    injectAperturePhot=True
    do_sim_posns=None#os.path.join(dataset,positionsf) # or None to not do this
    if useInjectedMap:
        perGalRadioMap=os.path.join(injectedMapDir,injectedMap)
    else:
        perGalRadioMap='/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS'
    perGalNoiseMap='/Users/jtlz2/video/vla/backrms.fits'

    if perGalUpsampleRate is None:
        upsampledMap=perGalRadioMap; upsampledNoise=perGalNoiseMap
    if perGalUpsampleRate is not None and perGalUpsampleRate < 0:
        if useInjectedMap:
            upsampledMap=os.path.join(injectedMapDir,injectedMap)
        else:
            upsampledMap='/Users/jtlz2/video/vla/VLA-VVDS_x%i.fits' % abs(perGalUpsampleRate)
    else:
        upsampledMap=None
    if perGalUpsampleRate is not None:
        upsampledNoise='/Users/jtlz2/video/vla/backrms_x%i.fits'% abs(perGalUpsampleRate)
    else:
        upsampledNoise=perGalNoiseMap
    perGalzScramble=False
    perGalzScrambleDZ=0.13
    perGalzScrambleDZSeed=1234
    perGalMScramble=False
    perGalMScrambleDM=0.1
    perGalMScrambleDMSeed=1235
    perGalXYScramble=False
    perGalXYRandomize=False#True # Randomize the position coordinates
    perGalXYRandomizeSeed=4321
    perGalXYRandomizeMask=False
    perGalXYRandomizeMaskArea=1.0#-(76123985.0/400320064.0)
    perGalXYScrambleN=10
    perGalXYScrambleM=40000
    perGalIDRange=[0,1e99]
    #perGalzBins=[-1.0e10,1.0e10]
    perGalzBins=[4.0,1.0e10]
    perGalMBins=[-100.0,100.0]
    perGalkabsBins=[-100.0,100.0]
    perGalNoiseThresh=1.0e10
    perGalClipThresh=1.0e10
    perGalCutKFlags=True
    perGalCutHaloFlags=True
    perGalKLim=23.5
    perGalzReadThresh=10.0 # Retain less than this
    perGalChiBestThresh=2.0e10
    perGalReliabilityThresh=-100.0 # Retain greater than this
    perGalCatForm=1
    perGalMasses='bon'
    perGalCat='pixels_Mz/pixels_%s_z_%3.1f_%3.1f_Ms_%3.1f_%3.1f_Kabs_%3.1f_%3.1f_noiseclip_1234.dat' % (perGalGalaxyType,perGalzBins[0],perGalzBins[-1],perGalMBins[0],perGalMBins[-1],perGalkabsBins[0],perGalkabsBins[-1])
    noise_only=None # uJy SURVEY_NOISE; None => use real map

    # Kim's cuts
    mcacuts=False
    if mcacuts:
        perGalIDRange=[0,2000000]
        perGalzBins=[-1.0e10,1.0e10]
        perGalMBins=[-100.0,100.0]
        perGalkabsBins=[-100.0,100.0]
        perGalNoiseThresh=1.0e10
        perGalClipThresh=1.0e10
        perGalCutKFlags=False
        perGalCutHaloFlags=False
        perGalKLim=23.4
        perGalzReadThresh=10.0 # Retain less than this
        perGalChiBestThresh=2.0e10
        perGalReliabilityThresh=0.8 # Retain greater than this
        perGalWrapRoot='video/%s' % perGalGalaxyType
        perGalCat='pixels_Mz/pixels_%s_z_%3.1f_%3.1f_Ms_%3.1f_%3.1f_Kabs_%3.1f_%3.1f_noiseclip_1234.dat' % (perGalGalaxyType,perGalzBins[0],perGalzBins[-1],perGalMBins[0],perGalMBins[-1],perGalkabsBins[0],perGalkabsBins[-1])

#-------------------------------------------------------------------------------

# Simulation parameters
SEED_SIM=1234

# Run simulate.py in batch mode
batch_sim=False
nbatch=1000

NLAWS_SIM=0
SLOPE_SIM=-2.0
##N_SIM=621114
N_SIM=None
##N_SIM=100000
C_SIM=26.0
SMIN_SIM=110.0
SMAX_SIM=839.3
NOISE_SIM=20.0
AREA_SIM=1.0
DUMP='fluxes.txt'
#DUMP=None

# FIRST large-area basis
#SLOPE_SIM=-1.5
#C_SIM=20.0
#SMIN_SIM=104.5
#SMAX_SIM=715.0
#NOISE_SIM=150.0
#AREA_SIM=100.0

# VIDEO basis
SLOPE_SIM=-1.94
C_SIM=10.0
SMIN_SIM=0.491
SMAX_SIM=853.0
NOISE_SIM=16.2
AREA_SIM=1.0

SLOPE_SIM=-1.95
C_SIM=10.0
SMIN_SIM=0.49
SMAX_SIM=850
NOISE_SIM=16.2
AREA_SIM=1.0

# VIDEO-style simulation (pre-flight run)
SLOPE_SIM=-1.95
#C_SIM=10.0
C_SIM=10.0
SMIN_SIM=0.49
SMAX_SIM=85.0
NOISE_SIM=16.2
AREA_SIM=1.0


# SKADS noise tests
#NOISE_SIM=0.5
SMIN_SKADS=0.01 # uJy
#SMIN_SKADS=0.5 # uJy
#SMASK_SKADS
SMAX_SKADS=85.0 # uJy
#SMAX_SKADS=600000.0 # uJy

SKADS_GO_VIA_MAP=True
NSKADS=None#72000 # or None to use all available sources for simulation
#NSKADS_RESCALING=373936.0/71962.0
#NSKADS_RESCALING=373924.0/71962.0
NSKADS_RESCALING=1.0

# VIDEO > 5-sigma style sim (pre-flight run)
#SLOPE_SIM=-1.95
#C_SIM=10.0
#SMIN_SIM=50.0
#SMAX_SIM=850.0
#NOISE_SIM=16.2
#AREA_SIM=1.0

#SMIN_SIM=0.49
#SMAX_SIM=1000.0


# VVDF style sim
#SLOPE_SIM=-2.0
#C_SIM=35.0
#SMIN_SIM=80.0
#SMAX_SIM=800.0
#NOISE_SIM=17.0
#AREA_SIM=1.0


# Ketron's Table 1
#C_SIM=40.0
#SLOPE_SIM=-1.50
#SMIN_SIM=1.00
#SMAX_SIM=20.00
#NOISE_SIM=10.0

D_SIM=-99.0
BETA_SIM=S0_SIM=GAMMA_SIM=S1_SIM=DELTA_SIM=S2_SIM=-99.0
#if NLAWS_SIM==0:
#    print 'Running SKADS sim'
if NLAWS_SIM==1:
    BETA_SIM=-99.0
    S0_SIM=-99.0
if NLAWS_SIM>1:
    #D_SIM=30.0
    BETA_SIM=-1.0
    S0_SIM=390.0
    S0_SIM=8.2
if NLAWS_SIM>2:
    GAMMA_SIM=-0.5
    S1_SIM=17.0
if NLAWS_SIM>3:
    DELTA_SIM=-0.3
    S2_SIM=40.0
    #print 'NLAWS > 2 not implemented for sims yet'

#if SIM_DO_CAT_NOISE:
#    OUTPUT='sim.txt'
#else:
#    OUTPUT='sim_noiseless.txt'

#-------------------------------------------------------------------------------

# Specify the data file (within outdir) and set up some survey parameters
numRedshiftSlices=len(whichRedshiftSlices)
numNoiseZones=len(noiseZones)
if dataset=='cosmos':
    #datafile='bondiover216_corr.txt'
    datafile='bondi2008_orig.txt'
    SURVEY_AREA=2.16 # sq. deg.
    SURVEY_NOISE=12.0 # uJy
elif dataset=='vvdf':
    datafile='bondi2003_orig.txt'
    SURVEY_AREA=1.00 # sq. deg.
    SURVEY_NOISE=17.0 # uJy
elif dataset=='video':
    #datafile='zwart2014_orig.txt'
    if binsHigh:
        datafile='%s_high.txt' % perGalGalaxyType
    else:
        datafile='%s_test_%i_%s.txt' % (perGalGalaxyType,nb,run_num_run)
    HALO_MASK=11436315.0/(19354.0*19354.0)
    SURVEY_AREA=1.0 *(1.0-HALO_MASK)# sq.deg. [Boris -> 0.97 sq. deg.]
    if perGalXYRandomizeMask: SURVEY_AREA=1.0
    SURVEY_NOISE=16.2 # uJy [median; mean=16.3]
elif 'sim' in dataset:
    if SKADS_GO_VIA_MAP:
        datafile='sim_extracted.txt'
    else:
        datafile='sim.txt'
    #datafile='sim_extracted_aper_inj.txt'
    #datafile='bondisim.txt'
    SURVEY_AREA=AREA_SIM # sq. deg.
    SURVEY_NOISE=NOISE_SIM # uJy
elif dataset == 'first':
    datafile='ketron2013_jz.txt'
    SURVEY_AREA=2.16 # sq. deg.
    SURVEY_NOISE=150.0 # uJy
elif dataset == 'mca':
    datafile='bondi2003_mca.txt'
    SURVEY_AREA=1.00
    SURVEY_NOISE=17.0
elif dataset == 'en1jvla':
    if not doPoln:
        stokes='I'
    else:
        stokes='P'
    weightvlaf='/Users/jtlz2/Dropbox/elaisn1/maps/EN1.I.mosaic.sensitivity_vla.fits'
    noisezonesf='noisezones_%s.txt'%stokes
    noisezonesf=os.path.join(dataset,noisezonesf)
    WEIGHT_COL=47
    cutsDict={'star':[30,0],'lacy':[34,-1],'stern':[38,-1],'donley':[40,-1],\
              'noise0':[WEIGHT_COL,1.0,1.4],'noise1':[WEIGHT_COL,1.4,1.96],\
              'noise2':[WEIGHT_COL,1.96,2.74],\
              'noise3':[WEIGHT_COL,2.74,3.84],'noise4':[WEIGHT_COL,3.84,5.38],\
              'noise5':[WEIGHT_COL,5.38,7.53],'noise6':[WEIGHT_COL,7.53,13.0]}
              #'noise0':[47,1.5,2.0]}
              #'noise0':[47,1.1,1.5]}

#    numNoiseZones=len([k for k in cutsDict.keys() if 'noise' in k])
    datafiles=['en1jvla_%s_a%i.txt'%(stokes,n) for n in noiseZones]
    datafile=datafiles[0]
    SURVEY_AREAS={}; SURVEY_NOISES={}
    if not os.path.exists(noisezonesf):
        nz=numpy.zeros((numNoiseZones,4))
    else:
        nz=numpy.atleast_2d(numpy.genfromtxt(noisezonesf))
    for i,d in enumerate(datafiles):
        #print nz
        SURVEY_AREAS[noiseZones[i]]=nz[noiseZones[i],-1]
        SURVEY_NOISES[noiseZones[i]]=nz[noiseZones[i],:-1]
    if not doPoln:
        SURVEY_NOISE=1.164 # uJy
    else:
        SURVEY_NOISE=1.034#0.8 # uJy
    SURVEY_AREA=sum([a for (f,a) in SURVEY_AREAS.iteritems()])
    redshifts=[0.0 for i in range(numNoiseZones)]
    #SURVEY_AREA=0.1
elif dataset=='sdss':
    SURVEY_NOISE=150.0 # uJy
    SURVEY_AREA=8609.67 # sq. deg.
    zmanifestf='sdss/sdss_manifest.txt'
    zmanifest=numpy.genfromtxt(zmanifestf)
    datafiles=['sdss_dr12s%i.txt'%i for i in zmanifest[[s-1 for s in whichRedshiftSlices],0]]
    datafile=datafiles[0]
    #numRedshiftSlices=zmanifest.shape[0]
    redshifts=zmanifest[[s-1 for s in whichRedshiftSlices],3]
    numRedshiftSlices=len(whichRedshiftSlices)

#-------------------------------------------------------------------------------

# Parameters for binning a catalogue
if context=='b' or True:
    BIN_CAT_CLIP=None
    CORR_RESOLUTION=None
    CORR_BINS=None
    BOUT_HISTO=os.path.join(dataset,'flux_histo_%s_%s.pdf'%(perGalGalaxyType,run_num))
    if dataset == 'first':
        #BIN_CAT='cats/5sigma.txt'
        BIN_CAT_FORM=0 # COSMOS-VLA
        BIN_CAT='cats/smap_scat.txt'
        BIN_COL_CLIP=1 # Retain only VLA-COSMOS CAT sources for which S >BIN_CAT_CLIP mJy
        #BIN_CAT_CLIP = 5.0 * 12.0 / 1000.0 # 5 sigma in mJy
        BIN_COL=0 # map-extracted [FIRST MAP]
        #BIN_COL=1 # cat-extracted [COSMOS 1.4-GHz CAT]
        BOUT_CAT='first/ketron2013_jz.txt'

    elif dataset == 'video':
        BIN_CAT_FORM=1 # VIDEO-VLA
        BIN_CAT='video/%s/%s'% (perGalGalaxyType,perGalCat)
        BIN_CAT_CLIP=None
        CORR_BINS=1.0/numpy.array([0.333,0.872,0.901,0.917,0.902,0.891,0.903,\
            0.889,1.000,0.875,0.667,1.000,0.500,1.000])
        CORR_BINS=1.0/(0.897*numpy.ones(nb)) # 41 or 24
        CORR_BINS=numpy.ones(nb) # 41 or 24
        if perGalXYRandomize:
            CORR_BINS=numpy.ones(nb)
        BIN_COL=12
        if binsHigh:
            BOUT_CAT='video/%s_high.txt' % perGalGalaxyType
        else:
            BOUT_CAT='video/%s_test_%i_%s.txt' % (perGalGalaxyType,nb,run_num)
        CORR_RESOLUTION=1.0/(1.0-0.024083)
        #CORR_RESOLUTION=1.20
        CORR_RESOLUTION=1.00
        CORR_RESOLUTION=1.00-0.92e-2

    elif dataset == 'vvdf':
        BIN_CAT_FORM=2 # VVDF format (web -> .xml -> .asc)
        BIN_CAT='vvdf/bondi2003.asc'
        BIN_CAT_CLIP=None
        BIN_COL=8 # [sic] - integrated flux, ignoring text COLUMNS
        BOUT_CAT='vvdf/bondi2003_jz.txt'
        CORR_BINS=numpy.array([1.29,1.25,1.00,1.00,1.00])

    elif dataset == 'mca':
        BIN_CAT_FORM=3
        BIN_CAT='mca/Likelihoodratio2013.csv'
        BIN_COL_CLIP=15
        BIN_CAT_CLIP=None # Cut on reliability
        CORR_BINS=1.0/numpy.array([0.333,0.872,0.901,0.917,0.902,0.891,0.903,\
                    0.889,1.000,0.875,0.667,1.000,0.500,1.000])
        CORR_BINS=1.0/(0.897*numpy.ones(14))
        BIN_COL=7 # [sic] - peak flux, ignoring text COLUMNS
        BOUT_CAT='mca/bondi2003_mca.txt'

    elif 'sim' in dataset:
        # Run directly from the injection_phot catalogue
        BIN_CAT_FORM=4
        #BIN_CAT='video/%s/%s'% (perGalGalaxyType,perGalCat)
        BIN_CAT_CLIP=None # Cut on reliability
        #CORR_BINS=numpy.ones((14))
        BIN_COL=12 # [sic] - peak flux, ignoring text COLUMNS
        BOUT_CAT=os.path.join(dataset,'sim_extracted.txt')
        CORR_RESOLUTION=1.0
        BIN_CAT=os.path.join(dataset,'injection_phot.txt')
        BIN_COL=7 #9 # [sic]

    elif dataset == 'en1jvla':
        BIN_CAT_FORM=8 # See /Users/jtlz2/elaisn1/servs/scripts/*.stil
        BIN_CAT='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-vla-overlap-with-detns-150813.txt'
        #BIN_CAT='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-vla-overlap-with-detns-shifted10p10p-150817.txt'
        #BIN_CAT='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-vla-overlap-with-detns-random-150817.txt'
        if doPoln:
            shiftedPolnPosns=True
            randomPolnPosns=False
            if shiftedPolnPosns:
                BIN_CAT='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-vla-overlap-with-detns-shifted10p10p-w-poln-151001.txt'
            elif randomPolnPosns:
                BIN_CAT='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-vla-overlap-with-detns-random-w-poln-151001.txt'
            else:
                BIN_CAT='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-vla-overlap-w-poln-with-detns-150818.txt'

        BIN_CAT_CLIP=None
        if doPoln:
            BIN_COL=49 # peak P flux
        else:
            BIN_COL=46 # peak I flux
        BOUT_CAT='en1jvla/en1jvla_%s.txt'%stokes

#-------------------------------------------------------------------------------

if context=='u' or context == 'i':
    upsamplingFactor=perGalUpsampleRate
    if upsamplingFactor is not None:
        upsamplingFactor=abs(perGalUpsampleRate)
        outputMap='/Users/jtlz2/video/vla/VLA-VVDS_x%i.fits' % upsamplingFactor
        outputNoiseMap='/Users/jtlz2/video/vla/backrms_x%i.fits' % upsamplingFactor
    injectInputMap='/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS'
    injectInputNoiseMap='/Users/jtlz2/video/vla/backrms.fits'
    if useInjectedMap and context != 'i':
        inputMap=os.path.join(dataset,injectedMap)
        outputMap=os.path.join(dataset,'sim_map_x2.fits')

if context=='i' or True:
    radioObservingFreqHz=1.4e9
    radioSynthBeamFWHM=4.0 # pixels/upsamplingFactor
    radioPixels2Arcsec=1.5
    radioSynthOmegaSr=sqDeg2sr*beamFac*(radioPixels2Arcsec*radioSynthBeamFWHM/3600.0)**2
    injectMaskRadius=0.5 # pixels/upsamplingFactor # is 0.5
    injectGaussianShiftX=0.0 # pixels
    injectGaussianShiftY=0.0 # pixels
    injectDoPosnNoise=False
    injectPositionNoise=0.1 # pixels/upsamplingFactor
    injectPhotometryDumpFile='injection_phot.txt'
    injectFakeNoise=True#False
    injectionNoise=NOISE_SIM*48.0240 # 48.0240 is for NOISE_SIM=16.2 uJy
    injectRecyclePositions=False
    injectPostageStampRadius=10.0 # pixels/upsamplingFactor
    injectUseSKADSPosns=True
    injectNumberOfSources=None # or None to inject all sources
    injectNumberOfSourcesExtracted=72000 # or None to extract all sources
    injectRandomizeExtractedPosns=True
    injectRandomizeExtractedPosnsSeed=4321
    injectDoBackSub=None#1.12 # None or flux to subtract in uJy
    injectSortFluxesDescending=True

#-------------------------------------------------------------------------------


# *** Avoid noiseless data
noisy=True

# STATIC: Fetch the data for this run (so they are available
#         independently of lumfunc; only if they're there)
# (These are the data for the realisation in Ketron's paper)
#ksRaw=numpy.array([122847,103812,87685,74436,63342,53420,45107,38055,32410])
#ksNoisy=numpy.array([8878,12240,17420,24516,34372,46473,60051,69770,63170])

# Need to sort this out
# Only really matters if simulating?
#datafile=os.path.join(outdir,datafile)
datafile=os.path.join(dataset,datafile)
for i,d in enumerate(datafiles):
    datafiles[i]=os.path.join(dataset,datafiles[i])

# This block is called nproc times - how to (circularly) use only the master?
# The thing to do is to read it from master and broadcast it
# Read in the data except in simulation mode
# The data are corrected by the completeness/ang. size factor and
#              normalized to 1 sq. deg.
if os.path.exists(datafile):
    if numNoiseZones==1 or numRedshiftSlices==0:
        # These are two hacks to change the recon binning scheme
        #if context=='r' and dataset=='first':
        #    datafile='%s/ketron2013_110_839_jz.txt' % dataset
        #if context=='r' and 'sim' in dataset:
        #    datafile='%s/sim_extracted.txt' % dataset
        data=numpy.genfromtxt(datafile)

        #if 'sim' in dataset or dataset in \
        #  ['cosmos','vvdf','first','video','mca','sdss','en1jvla']:
        bin_medians = data[:,2] # uJy [not that they are ever used -?]
        # In the table file, counts are for that SURVEY_AREA (so process THOSE)
        ksRaw       = data[:,3] * data[:,8] #/ (sqDeg2sr*SURVEY_AREA) #/ SURVEY_AREA # ???
        ksNoisy     = data[:,3] * data[:,8] #/ (sqDeg2sr*SURVEY_AREA) #/ SURVEY_AREA # ???

        # OLD VERSION SIM FORMAT:
        #elif dataset=='sim_old':
        #    bin_medians=data[:,0]; ksRaw=data[:,1]; ksNoisy=data[:,2]
            #print '***NB ksRaw -> junk'
            # Counts are now per sr

        # I retired these because a double import was applying them
        # twice (bug to be fixed...)
    #    print '-> Corrected counts for completeness and/or ang. size (C factor)'
        #print '-> units are sr^-1'
        #ksRaw *= data[:,8]
        #ksNoisy *= data[:,8]
        #print '-> Corrected counts for survey area (SURVEY_AREA)'
        #ksRaw *= 1.0/SURVEY_AREA
        #ksNoisy *= 1.0/SURVEY_AREA

        nRaw=numpy.sum(ksRaw)
        nNoisy=numpy.sum(ksNoisy)
        if noisy:
            ks=ksNoisy
        else:
            ks=ksRaw

else:
    if context=='y':
        print '***Cannot find an existing data file!! (%s)' % datafile
    if context in ['b','l']:
        x=os.path.join(dataset,'crash.txt')
        open(x, 'a').close()
        sys.exit(0)

#-------------------------------------------------------------------------------

# These are the truths
C_FIRST=7.2
SLOPE_FIRST=-2.11
SMIN_FIRST=104.5
SMAX_FIRST=715.0

C_RANDOM=4.5
SLOPE_RANDOM=-1.66
SMIN_RANDOM=1.03
SMAX_RANDOM=872.5

C_TARGET=19.7
SLOPE_TARGET=-2.32
SMIN_TARGET=110.3
SMAX_TARGET=839.3

SLOPE_VVDF=-2.28
CONVERT_VVDF=10**(-3.0*SLOPE_VVDF)*sqDeg2sr/1.0e3
C_VVDF=57.54/CONVERT_VVDF
SMIN_VVDF=80.0
SMAX_VVDF=600.0

SLOPE_VVDF2=-1.79
CONVERT_VVDF2=10**(-3.0*SLOPE_VVDF2)*sqDeg2sr/1.0e3
C_VVDF2=75.86/CONVERT_VVDF2
SMIN_VVDF2=600.0
SMAX_VVDF2=119230.0


#-------------------------------------------------------------------------------

# Set up the binning
#binstyle='sim'
#binstyle='kmw'
#binstyle='bondi2008'
#binstyle='bondi2003'
#binstyle='first'
#binstyle='video2014'
#binstyle='mca2014'
#binstyle='en1jvla'
binstyle='sdss'

#-------------------------------------------------------------------------------

# Set up the truths
C_TRUE=SLOPE_TRUE=SMIN_TRUE=SMAX_TRUE=BETA_TRUE=\
  S0_TRUE=GAMMA_TRUE=S1_TRUE=DELTA_TRUE=S2_TRUE=-99.0
if dataset == 'cosmos':
    C_TRUE=C_TARGET
    SLOPE_TRUE=SLOPE_TARGET
    SMIN_TRUE=SMIN_TARGET
    SMAX_TRUE=SMAX_TARGET
elif dataset == 'first':
    C_TRUE=C_FIRST
    SLOPE_TRUE=SLOPE_FIRST
    SMIN_TRUE=SMIN_FIRST
    SMAX_TRUE=SMAX_FIRST
elif dataset == 'vvdf':
    C_TRUE=C_VVDF
    SLOPE_TRUE=SLOPE_VVDF
    SMIN_TRUE=SMIN_VVDF
    SMAX_TRUE=SMAX_VVDF
elif 'sim' in dataset or 'kmw' in dataset:
    C_TRUE=C_SIM
    SLOPE_TRUE=SLOPE_SIM
    SMIN_TRUE=SMIN_SIM
    SMAX_TRUE=SMAX_SIM
    BETA_TRUE=BETA_SIM
    S0_TRUE=S0_SIM
    CONVERT_C_TRUE=10.0**(6.0*(SLOPE_SIM+2.5))
    C_TRUE=CONVERT_C_TRUE/C_SIM
    #print '-> Adjusted C_TRUE for SLOPE factor'

#-------------------------------------------------------------------------------


# Command-line arguments [mostly vestigial]
verbose=False
sigma=SURVEY_NOISE # uJy
area=SURVEY_AREA # sq. deg.
#model=1 # Vary all parameters
#loud=False
#tellthetruth=False


#-------------------------------------------------------------------------------

###REF_AREA=10.0 # sq. deg. for following arrays - vestigial - ***keep this fixed

# Hard-code the bins
if binstyle=='sim':
    bins=numpy.linspace(-270.0,1000.0,8)
    bins=numpy.array([-80.0,-50.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,20.0,\
                        50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])
    bins=numpy.array([-80.0,-50.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,20.0,\
                        50.0,80.0])
    #bins=numpy.array([-50.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,20.0,\
    #                    50.0,80.0])
    #bins=numpy.array([80.0,100.0,150.0,200.0,250.0,300.0,400.0,500.0,700.0,800.0])

    # This is to allow for different bins in the recon
    # Negative bins are now stripped out below
    ##if binsHigh or context=='r' :#or context=='s':
    #    bins=1000.0*numpy.array([0.1103,0.1351,\
    #                             0.1655,0.2028,0.2484,0.3043,0.3728,0.4566,\
    #                             0.5593,0.6851,0.8393]) # mJy -> uJy
    ##    bins=numpy.array([0.0,5.0,10.0,20.0,50.0,80.0])
    #bins=numpy.array([-100.0,-50.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,\
    #                  50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])
    #bins=numpy.array([5.0,10.0,20.0,\
    #                  50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])

    #bins=numpy.array([80.0,100.0,150.0,200.0,250.0,300.0,400.0,500.0,700.0,800.0])
    #bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #                         0.27,0.41,0.61]) # mJy -> uJy
    if nb==38:
        bins=numpy.array([-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,85.0])
        bins[0]=-69.0
    elif nb==37:
        bins=numpy.array([-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,85.0])
        bins[0]=-65.0
    elif nb==34:
        bins=numpy.array([-22.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0])
    elif nb==27:
        bins=numpy.array([-2.1,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0])
    elif nb==24:
        bins=numpy.array([-0.4,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0])
    elif nb==13:
        bins=numpy.array([-0.4,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0])

elif binstyle=='kmw':
    bins=numpy.array([1.000,1.394,1.945,2.714,3.786,5.281,7.368,10.27,14.33,20.00])
elif binstyle=='first':
    bins=numpy.linspace(-270.0,1000.0,8)
    # This is to allow for different bins in the recon
    if binsHigh or context=='r':
        bins=1000.0*numpy.array([0.1103,0.1351,\
                                 0.1655,0.2028,0.2484,0.3043,0.3728,0.4566,\
                                 0.5593,0.6851,0.8393]) # mJy -> uJy
    ##bins=1000.0*numpy.array([0.2028,0.2484,0.3043,0.3728,0.4566,\
    ##                         0.5593,0.6851]) # mJy -> uJy
    ##bins=numpy.linspace(-270.0,1000.0,8)
elif binstyle=='bondi2008':
    #bins=1000.0*numpy.array([0.0600,0.0735,0.0900,0.01103,0.1351,\
    #                         0.1655,0.2028,0.2484,0.3043,0.3728,0.4566,\
    #                         0.5593,0.6851,0.8393,1.0282]) # mJy -> uJy
    # 0.1103
    bins=1000.0*numpy.array([0.1103,0.1351,\
                             0.1655,0.2028,0.2484,0.3043,0.3728,0.4566,\
                             0.5593,0.6851,0.8393]) # mJy -> uJy
elif binstyle=='bondi2003':
    #bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #                         0.27,0.41,0.61,0.91,1.37,2.05,3.08,4.61,\
    #                         6.92,10.38,15.57,119.23]) # mJy -> uJy
    #bins=1000.0*numpy.array([0.01,0.12,0.18,\
    #                         0.27,0.41,0.80]) # mJy -> uJy
    bins=1000.0*numpy.array([0.08,0.12,0.18,\
                             0.27,0.41,0.61]) # mJy -> uJy

    #bins=1000.0*numpy.linspace(0.080,1.00,11)
    #bins=numpy.linspace(-470.0,839.0,7)
    #bins=numpy.linspace(80.0,610.0,5)
elif binstyle=='video2014' or binstyle=='mca2014':
    #bins=numpy.linspace(1.0,100.0,12) # uJy
    print '****Warning: Tailoring bins to galaxy type'
    if perGalGalaxyType in ['ell','sbn']:
        if nb==38:
            bins=numpy.array([-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
            bins[0]=-67.0
            if perGalGalaxyType=='sbn':
                bins[0]=-69.0
    if perGalGalaxyType in ['all','mid']:
        bins=numpy.linspace(-100.0,800.0,25) # uJy --- all
        bins=numpy.array([-79.0,-50.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,\
                          50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])
        bins=numpy.array([-108.0,-50.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,20.0,\
                          50.0,85.0])
        if nb==24:
            bins=numpy.array([-108.0,-80.0,-50.0,-30.0,-20.0,-10.0,-5.0,-2.5,-1.0,-0.5,0.0,0.5,1.0,2.5,5.0,7.5,10.0,15.0,20.0,25.0,30.0,40.0,50.0,65.0,85.0])
        elif nb==41 and not perGalXYRandomize:
            bins=numpy.array([-108.0,-90.0,-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
        elif nb==40 and not perGalXYRandomize:
            bins=numpy.array([-104.0,-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
        elif nb==39 and not perGalXYRandomize:
            bins=numpy.array([-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
            bins[0]=-73.0
        elif nb==38 and not perGalXYRandomize:
            bins=numpy.array([-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
            bins[0]=-57.0
        elif nb==37 and not perGalXYRandomize:
            bins=numpy.array([-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
            bins[0]=-57.0
        elif nb==31 and not perGalXYRandomize:
            bins=numpy.array([-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,85.0])
        elif nb==26 and not perGalXYRandomize:
            bins=numpy.array([-67.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,0.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,85.0])

        elif nb==50 and perGalXYRandomize:
            bins=numpy.array([-83.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0,110.3,135.1,165.5,202.8,248.4,304.3,372.8,456.6,559.3,685.1,839.3])
            #bins[0]=-116.0
        elif nb==51 and not perGalXYRandomize:
            bins=numpy.array([-81.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0,110.3,135.1,165.5,202.8,248.4,304.3,372.8,456.6,559.3,685.1,839.3])
        elif nb==52 and not perGalXYRandomize:
            bins=numpy.array([-108.0,-90.0,-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0,110.3,135.1,165.5,202.8,248.4,304.3,372.8,456.6,559.3,685.1,839.3])
        #bins[0]=-100.0
        elif nb==39 and perGalXYRandomize:
            bins=numpy.array([-67.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
            bins[0]=-83.0
        elif nb==40 and perGalXYRandomize:
            bins=numpy.array([-86.0,-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
            bins[0]=-91.0
        elif nb==41 and perGalXYRandomize:
            bins=numpy.array([-92.0,-90.0,-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])
        elif nb==55:
            bins=numpy.concatenate((numpy.array([-108.0,-90.0,-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0]),1000.0*numpy.array([0.12,0.18,0.27,0.41,0.61,0.91,1.37,2.05,3.08,4.61,6.92,10.38,15.57,119.23])))
        elif nb==37 and noise_only is not None:
            bins=numpy.array([-64.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0])

    if False and perGalGalaxyType=='ell':
        if nb==38:
            bins=numpy.array([-66.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])

    if False and perGalGalaxyType=='sbn':
        if nb==38:
            bins=numpy.array([-69.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,80.0,85.0])

elif binstyle=='en1jvla':
    #bins=numpy.linspace(-4.6,5.8,40)
    #bins=numpy.linspace(-15.0,5.8,40)
    #bins=numpy.linspace(-5.3,5.8,32)
    binsZoned={}

    # For shifted data
    if shiftedIntensityPosns:
        binsZoned[0]=numpy.array([-4.1,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[1]=numpy.array([-6.1,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[2]=numpy.array([-8.1,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.02,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])

    # For random-position data
    elif randomIntensityPosns:
        binsZoned[0]=numpy.array([-4.4,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[1]=numpy.array([-7.2,-6.5,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[2]=numpy.array([-7.2,-6.5,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])

    else:
        #binsZoned[0]=numpy.array([-5.3,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,-0.005,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        #binsZoned[1]=numpy.array([-9.5,-6.5,-5.3,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        #binsZoned[0]=numpy.array([-4.8,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        #binsZoned[1]=numpy.array([-6.3,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        #binsZoned[2]=numpy.array([-7.6,-6.5,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[0]=numpy.array([-4.7,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[1]=numpy.array([-5.7,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[2]=numpy.array([-6.7,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])

        binsZoned[0]=numpy.array([-4.8,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.005,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[1]=numpy.array([-6.3,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[2]=numpy.array([-8.9,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[3]=numpy.array([-12.3,-10.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[4]=numpy.array([-17.6,-15.0,-10.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[5]=numpy.array([-22.6,-20.0,-15.0,-10.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        binsZoned[6]=numpy.array([-36.3,-25.0,-20.0,-15.0,-12.0,-10.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.5,-0.1,-0.05,0.0,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
    if doPoln:
        #binsZoned[0]=numpy.array([0.0,0.125,0.25,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        #binsZoned[1]=numpy.array([0.0,0.125,0.25,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        #binsZoned[2]=numpy.array([0.0,0.125,0.25,0.5,1.0,2.0,3.0,4.0,5.0,5.8])
        if shiftedPolnPosns:
            binsZoned[0]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,5.8]) # 4.5
            binsZoned[1]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,5.8]) # 4.5

        elif randomPolnPosns:
            binsZoned[0]=numpy.array([0.03,0.05,0.08,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,5.8])
            binsZoned[0]=numpy.array([0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0])
            binsZoned[1]=numpy.array([0.03,0.05,0.08,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,5.8])
            if maskSERVS or maskSERVSplusPOLN:
                binsZoned[0]=numpy.array([0.015,0.03,0.05,0.08,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0])#5.8])
                
        else:
            binsZoned[0]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0])#,5.8])
            binsZoned[1]=numpy.array([0.03,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,5.8])
            binsZoned[2]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                                      0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,5.8])

        #binsZoned[1]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
        #                          0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0,5.5,5.8])
        #binsZoned[2]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
        #                          0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0,5.5,5.8])
        #binsZoned[1]=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
        #                          0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,4.0,5.0,5.5,5.8,6.5,7.0])

    bins=binsZoned[noiseZones[0]]

elif binstyle=='goodsnx':
    binsZoned={}
    binsZoned[0]=numpy.array([-4.2,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.75,-0.5,-0.3,-0.2,-0.1,-0.05,-0.02,-0.01,0.0,0.01,0.02,0.05,0.1,0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5])
    if doPoln:
        binsZoned[0]=numpy.array([0.0,0.02,0.05,0.1,0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5])
        # For shifted or random data
        #binsZoned[0]=numpy.array([0.0,0.03,0.05,0.1,0.3,0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5])

    assert(len(noiseZones)==1), 'too many noise zones'
    bins=binsZoned[noiseZones[0]]

elif binstyle=='sdss':
    bins=[ -700., -340,  -295,  -251,  -206,\
        -161. ,-138., -116., -93. ,-71., -48.,\
        -26. ,  -4. ,18., 40.,   63., 85.,  108.,\
         130,  153, 175,  197., 219 , 242., 264.,\
         287., 309,  332., 354., 377., 422., 467.,\
         512., 557., 646.,736.,   826.,  1006.,\
         1140. , 1320.,  1500.]
    
    #bins=bins[:38]
    #assert(len(bins)-1==nb)

    #bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #                         0.27,0.41,0.61]) # mJy -> uJy

    #bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #                         0.27,0.41,0.61,0.91,1.37,2.05,3.08,4.61,\
    #                         6.92,10.38,15.57,119.23]) # mJy -> uJy


    #elif perGalGalaxyType=='ell':
    #    bins=numpy.array([-79.0,-50.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,\
    #                      50.0,100.0,200.0,300.0,400.0,800.0])
    #elif perGalGalaxyType=='mid':
    #    bins=numpy.linspace(-100.0,800.0,25) # uJy --- mid
    #    bins=numpy.array([-79.0,-50.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,\
    #                      50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])
    #elif perGalGalaxyType=='sbn':
    #    bins=numpy.linspace(-100.0,800.0,25) # uJy --- mid
    #    bins=numpy.array([-79.0,-50.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,\
    #                      50.0,100.0,200.0,300.0,800.0])

    #bins=numpy.linspace(-100.0,450.0,12) # uJy --- ell, mid, sbn
    #if binsHigh or context=='r':
    #    bins=numpy.logspace(1.0,3.0,10) # 10**1.0 -> 10**3.0
    #    if context=='r':
    #        bins=numpy.logspace(numpy.log10(SURVEY_NOISE),3.0,10)
    #    if context=='b':
    #        bins=numpy.logspace(numpy.log10(5.0*SURVEY_NOISE),3.0,10)
    #    bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #                            0.27,0.41,0.61]) # mJy -> uJy
    #    bins=1000.0*numpy.array([0.08,0.12,0.18,0.27,0.41,0.61,\
    #                             0.91,1.37,2.05,3.08,4.61,6.92,10.38,15.57,119.23])
    #    bins=numpy.array([0.0,5.0,10.0,20.0,50.0,85.0])
    #    bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #                            0.27,0.41,0.61,0.91]) # mJy -> uJy
    #    #bins=1000.0*numpy.array([0.08,0.12,0.18,\
    #    #                        0.27,0.41,0.61,0.91,1.37,2.05,3.08,4.61,\
    #    #                        6.92,10.38,15.57,119.23]) # mJy -> uJy

# Eliminate negative bins for the reconstruction plot
if binsHigh or context=='r':
    bins_saved=bins
    bins=bins[numpy.where(bins>=0.0)]

nbins=len(bins)
dbins=numpy.gradient(bins)
# (bin medians now below the medianArray function, below)


#-------------------------------------------------------------------------------

# Priors

#C_MIN=1.0
#C_MAX=300.0
#SLOPE_MIN=-2.5
#SLOPE_MAX=-0.1
#SMIN_MIN=1.0
#SMIN_MAX=500.0
#SMAX_MIN=500.0
#SMAX_MAX=1000.0
##SMAX_MIN=10.0
##SMAX_MAX=50.0

#C_PRIOR='U'
C_PRIOR='LOG'

C_MIN=1.0             # sr^-1
if 'first' in dataset:
    C_MAX=1.0e5      # FIRST prior
    #C_MAX=10.0      # FIRST prior
elif nlaws==2:
    C_MAX=100.0       # TPL prior
elif binstyle=='kmw':
    C_MAX=1.0e6       # Ketron Table 1 prior
else:
    C_MAX=5000.0
C_MAX=100000.0

SLOPE_MIN=-2.5
SLOPE_MAX=-0.1

#C_MIN=C_MAX=C_TRUE
#SLOPE_MIN=SLOPE_MAX=SLOPE_TRUE

#SMIN_MIN=1.0        # uJy
#SMIN_MAX=300.0
#SMAX_MIN=400.0      # uJy
#SMAX_MAX=1000.0

#SMIN_MIN=1.0        # uJy
#SMIN_MAX=50.0
#SMAX_MIN=50.0      # uJy
#SMAX_MAX=150.0

# FIRST, COSMOS priors
SMIN_MIN=1.0        # uJy
SMIN_MAX=300.0
SMAX_MIN=500.0      # uJy
SMAX_MAX=1000.0

# FIRST priors
C_PRIOR='U'
C_MIN=1.0e-5             # sr^-1
C_MAX=100.0
SMIN_MIN=1.0        # uJy
SMIN_MAX=500.0
SMAX_MIN=500.0      # uJy
SMAX_MAX=1000.0

#C_PRIOR='U'
#C_MAX=1.0e6
#SMIN_MIN=SMIN_MAX=SMIN_TRUE       # delta fn on Smin
#SMAX_MIN=SMAX_MAX=SMAX_TRUE      # delta fn on Smax
#SLOPE_MIN=SLOPE_MAX=SLOPE_TRUE
#C_MIN=C_MAX=C_TRUE

# Ketron Table 1 priors
#SMIN_MIN=0.01        # uJy
#SMIN_MAX=5.0
#SMAX_MIN=10.0      # uJy
#SMAX_MAX=50.0


# VIDEO 2014 priors
# PRIORS_xyz
##SMIN_MIN=0.01        # uJy
##SMIN_MAX=5.0
##SMAX_MIN=16.0      # uJy
##SMAX_MAX=80.0
##C_PRIOR='LOG'
##C_MIN=1.0
##C_MAX=1000.0
#SMIN_MIN=0.01        # uJy
#SMIN_MAX=20.0 # Was 20.0
##SMIN_MIN=0.05        # uJy
##SMIN_MAX=1.0
#SMAX_MIN=20.0      # uJy
#SMAX_MAX=100.0

# ELAIS-N1 JVLA priors
SMIN_MIN=0.001 # i.e. 1 nJy
SMIN_MAX=1.5
SMAX_MIN=1.5
SMAX_MAX=5.8

# Based on S0=30 uJy - give SMAX the freedom to be roughly there
#SMAX_MIN=20.0      # uJy
#SMAX_MAX=30.0

# SKADS noise tests
#SMIN_MIN=0.01        # uJy
#SMIN_MAX=1.0
#SMAX_MIN=1.0      # uJy
#SMAX_MAX=100.0

#SMAX_MAX=200.0
#SMAX_MAX=100.0

C_PRIOR='LOG'
C_MIN=1.0e-5
C_MAX=1.0e9

# VIDEO > 5 sigma priors
#SMIN_MIN=0.01        # uJy
#SMIN_MAX=100.0
#SMAX_MIN=100.0      # uJy
##SMAX_MAX=1000.0
#SMAX_MAX=1000.0
#C_PRIOR='LOG'
#C_MIN=1.0e-5
#C_MAX=1.0e3


#if dataset=='video' and perGalGalaxyType=='ell':
#    SMIN_MAX=5.0
#if dataset=='video' and perGalGalaxyType=='sbn':
#    C_MAX=100.0
#if dataset=='video' and binsHigh:
#    C_MAX=100.0
#    SMIN_MIN=5.0*SURVEY_NOISE        # uJy
#    SMIN_MAX=500.0
#    SMAX_MIN=500.0      # uJy
#    SMAX_MAX=1000.0

# VVDF priors
#SMIN_MIN=1.0        # uJy
##SMIN_MAX=500.0
#SMIN_MAX=300.0
#SMAX_MIN=500.0      # uJy
#SMAX_MAX=1000.0
#C_PRIOR='LOG'
#C_MIN=1.0e-5
##C_MAX=1000.0
#C_MAX=500.0



#SMIN_MIN=1.0       # uJy
#SMIN_MAX=40.0
#SMAX_MIN=40.0        # uJy
#SMAX_MAX=100.0


#SMIN_MIN=SMIN_MAX=1.0       # delta fn on Smin
#SMAX_MIN=SMAX_MAX=20.0      # delta fn on Smax
#C_MIN=C_MAX=40.0            # delta fn on C
#SLOPE_MIN=SLOPE_MAX=-1.5    # delta fn on slope

#SMIN_MIN=SMIN_MAX=110.0       # delta fn on Smin
#SMAX_MIN=SMAX_MAX=839.0      # delta fn on Smax
#C_MIN=C_MAX=19.7            # delta fn on C
#SLOPE_MIN=SLOPE_MAX=-2.32    # delta fn on slope

#======================================================
# delta-fn settings:
#SMIN_MIN=SMIN_MAX=SMIN_TRUE       # delta fn on Smin
#SMAX_MIN=SMAX_MAX=SMAX_TRUE      # delta fn on Smax
#C_MIN=C_MAX=C_TRUE            # delta fn on C
#SLOPE_MIN=SLOPE_MAX=SLOPE_TRUE    # delta fn on slope
#======================================================

# Priors for polynomial coefficients
POLYCOEFF_PRIOR='LOG'
#POLYCOEFF_MIN=1.0e9
#POLYCOEFF_MAX=1.0e23
POLYCOEFF_MIN=1.0e5
POLYCOEFF_MAX=1.0e24

# Priors for bin/pole/node amplitudes
POLEAMPS_PRIOR='LOG'
POLEAMPS_MIN=1.0e3
POLEAMPS_MAX=1.0e10


#D_MIN=1.0
#D_MAX=100.0
BETA_MIN=-2.5
BETA_MAX=-0.1
S0_MIN=SMIN_MIN
S0_MAX=SMAX_MAX

#S0_MIN=SMIN_MAX
#S0_MAX=SMAX_MIN
GAMMA_MIN=-2.5
GAMMA_MAX=-0.1
S1_MIN=SMIN_MIN
S1_MAX=SMAX_MAX

DELTA_MIN=-2.5
DELTA_MAX=-0.1
S2_MIN=SMIN_MIN
S2_MAX=SMAX_MAX

#S0_MIN=S0_MAX=15.0
#SMIN_MIN=SMIN_MAX=1.21
#SMIN_MIN=1.21


# LF parameters

# Prior types
LMIN_PRIOR='LOG'
LMAX_PRIOR='LOG'
LSTAR_PRIOR='LOG'
LSLOPE_PRIOR='U'
LNORM_PRIOR='LOG'
LZEVOL_PRIOR='U'
LSLOPE2_PRIOR='U'

# Prior ranges
LMIN_MIN=1.0e20
LMIN_MAX=1.0e23
LMAX_MIN=1.0e23
LMAX_MAX=1.0e25
LSTAR_MIN=1.0e20
LSTAR_MAX=1.0e25
LSLOPE_MIN=0.0
LSLOPE_MAX=-8.0
LNORM_MIN=1.0e-10
LNORM_MAX=1.0e-5
LSLOPE2_MIN=0.0
LSLOPE2_MAX=-8.0

# LF redshift evolution (if > 1 z slice only)
LZEVOL_MIN=-5.0
LZEVOL_MAX=5.0
if len(whichRedshiftSlices) <= 1:
    LZEVOL_MIN=LZEVOL_MAX=0.0 # delta fn on LZEVOL


# MF parameters

# Prior types
MMIN_PRIOR='LOG'
MMAX_PRIOR='LOG'
MSTAR_PRIOR='LOG'
MSLOPE_PRIOR='U'
MNORM_PRIOR='LOG'
MZEVOL_PRIOR='U'
MSLOPE2_PRIOR='U'

# Prior ranges
MMIN_MIN=1.0e20
MMIN_MAX=1.0e23
MMAX_MIN=1.0e23
MMAX_MAX=1.0e25
MSTAR_MIN=1.0e20
MSTAR_MAX=1.0e25
MSLOPE_MIN=0.0
MSLOPE_MAX=-8.0
MNORM_MIN=1.0e-10
MNORM_MAX=1.0e-5
MSLOPE2_MIN=0.0
MSLOPE2_MAX=-8.0

# MF redshift evolution (if > 1 z slice only)
MZEVOL_MIN=-5.0
MZEVOL_MAX=5.0
if len(whichRedshiftSlices) <= 1:
    MZEVOL_MIN=MZEVOL_MAX=0.0 # delta fn on MZEVOL


#S0_MIN=S0_MAX=S0_TRUE            # delta fn on S0
#BETA_MIN=BETA_MAX=BETA_TRUE    # delta fn on beta
#D_MIN=D_MAX=D_TRUE            # delta fn on C
    
if numNoiseZones>1:
    NOISES_MIN={}; NOISES_MAX={}
    for nn in range(numNoiseZones):
        NOISES_MIN[noiseZones[nn]]=SURVEY_NOISES[noiseZones[nn]][1]
        NOISES_MAX[noiseZones[nn]]=SURVEY_NOISES[noiseZones[nn]][2]
else:
    #NOISE_MIN=0.5*SURVEY_NOISE
    #NOISE_MAX=2.0*SURVEY_NOISE
    NOISE_MIN=0.5*SURVEY_NOISES[noiseZones[0]][1]
    NOISE_MIN=SURVEY_NOISES[noiseZones[0]][1]
    #NOISE_MIN=1.0 # Fixed
    NOISE_MAX=2.0*SURVEY_NOISES[noiseZones[0]][2]
    NOISE_MAX=SURVEY_NOISES[noiseZones[0]][2]

#NOISE_MIN=NOISE_MAX=SURVEY_NOISE
if not floatNoise:
    #NOISE_MIN=NOISE_MAX=SURVEY_NOISE # delta fn on NOISE
    NOISE_MIN=NOISE_MAX=(SURVEY_NOISES[noiseZones[0]][1]+SURVEY_NOISES[noiseZones[0]][2])/2.0 # delta fn on NOISE

if doPoln and doRayleigh:
    #C_PRIOR='LOG'
    #C_MIN=C_MAX=0.0
    SMIN_MIN=0.001#1.0e-12
    SMIN_MAX=0.001
    SMAX_MIN=SMAX_MAX=1.5
    #SMAX_MIN=SMAX_MAX=binsZoned[noiseZones[0]][-1]#1.0e12
    S0_MIN=S0_MAX=2.0e-3
    S1_MIN=S1_MAX=3.0e-3
    S2_MIN=S2_MAX=4.0e-3
    SLOPE_MIN=SLOPE_MAX=-2.0
    #C_MIN=CMAX=1.0e-50

#-------------------------------------------------------------------------------

# Experiment resolution (for source-confusion prior/s)
confusionNoiseThreshold=0.5 # x sigma
checkConfusionNoise=False
confusionNoiseCheckSmax=5.0*SURVEY_NOISE # or None => Smax for each itern

checkSmin=False
numericalCheckSmin=True
dOmega=7.72e-6 # This is 10 x 10 arcsec^2 / sq. deg.
#dOmega=7.72e-7 # This is 10 arcsec^2 / sq. deg.
#dOmega=2.44e-7 
#dOmega=2.235e-10
#dOmega=2.7778e-3
# Confusion-noise factor
div=2.0
# Non-overlapping prior on Smin/max
assert(SMAX_MIN >= SMIN_MAX), 'Smin/Smax priors must not overlap!'
#if nlaws > 1:
#    assert(S0_MIN >= SMIN_MAX), 'Smin/S0 priors must not overlap!'
#    assert(SMAX_MIN >= S0_MAX), 'S0/Smax priors must not overlap!'

#-------------------------------------------------------------------------------

# Set up the parameters for triangle plots and reconstruction


if nlaws == 1:
    parameters=['C','slope','Smin','Smax']
    plotRanges={'C':[C_MIN,C_MAX],
                'slope':[SLOPE_MIN,SLOPE_MAX],
                'Smin':[SMIN_MIN,SMIN_MAX],
                'Smax':[SMAX_MIN,SMAX_MAX]}
    plotTruth={'C':C_TRUE,
               'slope':SLOPE_TRUE,
               'Smin':SMIN_TRUE,
               'Smax':SMAX_TRUE}
elif nlaws == 2:
    parameters=['C','slope','Smin','Smax','beta','S0']
    plotRanges={'C':[C_MIN,C_MAX],
                'slope':[SLOPE_MIN,SLOPE_MAX],
                'Smin':[SMIN_MIN,SMIN_MAX],
                'Smax':[SMAX_MIN,SMAX_MAX],
                'beta':[BETA_MIN,BETA_MAX],
                'S0':[S0_MIN,S0_MAX]}
    plotTruth={'C':C_TRUE,
               'slope':SLOPE_TRUE,
               'Smin':SMIN_TRUE,
               'Smax':SMAX_TRUE,
               'beta':BETA_TRUE,
               'S0':S0_TRUE}
elif nlaws == 3:
    parameters=['C','slope','Smin','Smax','beta','S0','gamma','S1']
    plotRanges={'C':[C_MIN,C_MAX],
                'slope':[SLOPE_MIN,SLOPE_MAX],
                'Smin':[SMIN_MIN,SMIN_MAX],
                'Smax':[SMAX_MIN,SMAX_MAX],
                'beta':[BETA_MIN,BETA_MAX],
                'S0':[S0_MIN,S0_MAX],
                'gamma':[GAMMA_MIN,GAMMA_MAX],
                'S1':[S1_MIN,S1_MAX]}
    plotTruth={'C':C_TRUE,
               'slope':SLOPE_TRUE,
               'Smin':SMIN_TRUE,
               'Smax':SMAX_TRUE,
               'beta':BETA_TRUE,
               'S0':S0_TRUE,
               'gamma':GAMMA_TRUE,
               'S1':S1_TRUE}
elif nlaws == 4:
    parameters=['C','slope','Smin','Smax','beta','S0','gamma','S1','delta','S2']
    plotRanges={'C':[C_MIN,C_MAX],
                'slope':[SLOPE_MIN,SLOPE_MAX],
                'Smin':[SMIN_MIN,SMIN_MAX],
                'Smax':[SMAX_MIN,SMAX_MAX],
                'beta':[BETA_MIN,BETA_MAX],
                'S0':[S0_MIN,S0_MAX],
                'gamma':[GAMMA_MIN,GAMMA_MAX],
                'S1':[S1_MIN,S1_MAX],
                'delta':[DELTA_MIN,DELTA_MAX],
                'S2':[S2_MIN,S2_MAX]}
    plotTruth={'C':C_TRUE,
               'slope':SLOPE_TRUE,
               'Smin':SMIN_TRUE,
               'Smax':SMAX_TRUE,
               'beta':BETA_TRUE,
               'S0':S0_TRUE,
               'gamma':GAMMA_TRUE,
               'S1':S1_TRUE,
               'delta':DELTA_TRUE,
               'S2':S2_TRUE}
elif modelFamily=='LFsch':
    parameters=['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL']
    plotRanges={'LNORM':[LNORM_MIN,LNORM_MAX],
                'LSTAR':[LSTAR_MIN,LSTAR_MAX],
                'LSLOPE':[LSLOPE_MIN,LSLOPE_MAX],
                'LMIN':[LMIN_MIN,LMIN_MAX],
                'LMAX':[LMAX_MIN,LMAX_MAX],
                'LZEVOL':[LZEVOL_MIN,LZEVOL_MAX]}
    plotTruth={'LNORM':LNORM_TRUE,
               'LSTAR':LSTAR_TRUE,
               'LSLOPE':LSLOPE_TRUE,
               'LMIN':LMIN_TRUE,
               'LMAX':LMAX_TRUE,
               'LZEVOL':LZEVOL_TRUE}
elif modelFamily=='LFdpl':
    parameters=['LNORM','LSTAR','LSLOPE','LSLOPE2','LMIN','LMAX','LZEVOL']
    plotRanges={'LNORM':[LNORM_MIN,LNORM_MAX],
                'LSTAR':[LSTAR_MIN,LSTAR_MAX],
                'LSLOPE':[LSLOPE_MIN,LSLOPE_MAX],
                'LSLOPE':[LSLOPE2_MIN,LSLOPE2_MAX],
                'LMIN':[LMIN_MIN,LMIN_MAX],
                'LMAX':[LMAX_MIN,LMAX_MAX],
                'LZEVOL':[LZEVOL_MIN,LZEVOL_MAX]}
    plotTruth={'LNORM':LNORM_TRUE,
               'LSTAR':LSTAR_TRUE,
               'LSLOPE':LSLOPE_TRUE,
               'LSLOPE2':LSLOPE2_TRUE,
               'LMIN':LMIN_TRUE,
               'LMAX':LMAX_TRUE,
               'LZEVOL':LZEVOL_TRUE}
elif modelFamily=='poly':
    parameters=['p%i' % ic for ic in range(nlaws)]
    for ic in range(nlaws):
        plotRanges['p%i'%ic] = [POLYCOEFF_MIN,POLYCOEFF_MAX]
        plotTruth['p%i'%ic] = -1000.0

if floatNoise:
    if numNoiseZones>1:
        for nn in range(numNoiseZones):
            parameters += ['noise%i'%nn]
            plotRanges['noise%i'%nn]=[NOISES_MIN[noiseZones[nn]],\
                                      NOISES_MAX[noiseZones[nn]]]
            plotTruth['noise%i'%nn]=SURVEY_NOISE
    else:
        parameters += ['noise']
        plotRanges['noise']=[NOISE_MIN,NOISE_MAX]
        plotTruth['noise']=SURVEY_NOISE


# Triangle plot
triangle='triangle_%s.png' % run_num

# Plotting options - reconstruction
POINTS_OFFSET=0.0
#PLOT_XMIN=1.0      # uJy
#PLOT_XMAX=10000.0  # uJy
PLOT_LABEL=''
labelDict={'C':r'$C/$Jy$^{-1}$sr$^{-1}$',\
           'slope':r'$\alpha$','Smin':r'$S_{\mathrm{min}}$/$\mu$Jy',\
           'Smax':r'$S_{\mathrm{max}}$/$\mu$Jy','beta':r'$\beta$',\
           'S0':r'$S_0/\mu$Jy','gamma':r'$\gamma$','S1':r'$S_1/\mu$Jy',\
           'delta':r'$\delta$','S2':r'$S_2/\mu$Jy','sigma':r'$\sigma/\mu$Jy',\
           'LNORM':r'$\rho^*$','LSTAR':r'$L^*/WHz^{-1}$','LSLOPE':r'$\alpha_{LF}$',\
           'LMIN':r'$L_{\mathrm{min}}/WHz^{-1}$','LMAX':r'$L_{\mathrm{max}}/WHz^{-1}$/',\
           'LZEVOL':r'$\beta_{LF}$','LSLOPE2':r'$\alpha_{LF,2}$'}

#-------------------------------------------------------------------------------

# Set some MultiNEST parameters here
outstem='1-'

n_live_points=500 # was 1000
multimodal=False
max_modes=3
SEED_SAMP=SEED_SAMP # [-1 for clock]


# Switch to INS
do_INS=False
#n_live_points=500

max_iter=0
evidence_tolerance=0.5 # was 0.5
sampling_efficiency=0.8 # default is 0.8

# Warning messages etc.
#print 'MOTD: %s' % MOTD

#laws={1:'SPL',2:'TPL',3:'XPL',4:'QPL',0:'SKADS'}
#print '****Considering %s model applied to %s law' % (laws[nlaws],laws[NLAWS_SIM])

#-------------------------------------------------------------------------------

#NOISE_SIM=None
