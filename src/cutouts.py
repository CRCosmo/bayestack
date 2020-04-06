from pylab import*
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
import numpy,random
from utils import gaussian
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Circle, Rectangle

data = fits.getdata('/home/eliab/Documents/OC/Project/Cosmos/vla_3ghz_msmf.fits')
h = fits.getheader('/home/eliab/Documents/OC/Project/Cosmos/vla_3ghz_msmf.fits')
w = wcs.WCS(h)

data2= data.reshape(data.shape[2:])
print 'reshaped'
ra,dec,sep = numpy.loadtxt('cosmos_d6_0_8arc.txt',unpack=True, usecols=[0,1,-1])

n=len(ra)/10
print n,len(ra)
size=2000

for i in range(1):
    index =random.randint(1,500)
    print 'start'
    hdr=None
    co=None

    print ra[index], dec[index]
    #continue

    r,d,o,p=w.wcs_world2pix(ra[index],dec[index],0,1.501191666680E+02,2.205833333330)
    print round(r),round(d)
    co = Cutout2D(data2,position=(round(r)-1, round(d) -1),size=size)
    print 'done cutting'
    
    print h['CRPIX1']
    print h['NAXIS1']
    print h['NAXIS2']
    print h['CRPIX2']
    
    
    print 'change'
    h['NAXIS']=4
   
    h['NAXIS1']=size
    h['NAXIS2']=size
    
    h['CRPIX1']=size/2
    h['CRPIX2']=size/2
    h['CRVAL1']=ra[index]
    h['CRVAL2']=dec[index]


    print h['CRPIX1']
    print h['NAXIS1']
    print h['NAXIS2']
    print h['CRPIX2']

    
    hdr=fits.PrimaryHDU(data=co.data,header=h)
    
    
    print 'fits created'  
    
    k=True
    s=0

    while k:
        try:
            hdr.writeto('%s.fits'%(s+1))
            k=False
        except IOError:
            print s
            s+=1            
        
    print 'fits saved'    
    
