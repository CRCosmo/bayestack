from pylab import*
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
import numpy,random,sys
from utils import gaussian
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Circle, Rectangle
import matplotlib.colors as colors


#data = fits.getdata('1.fits')
#h = fits.getheader('1.fits')
dfile='1.fits'
data = fits.getdata(dfile)
h = fits.getheader(dfile)

w = wcs.WCS(h)

#ra,dec,sep = numpy.loadtxt('cosmos_d2_0_8arc.txt',unpack=True, usecols=[0,1,-1])

print h['CRVAL1']

ax=subplot()

for i in range(9):
 ra,dec,sep = numpy.loadtxt('cos_data/cosmos_d%d_0_8arc.txt'%(i+1),unpack=True, usecols=[0,1,-1])
 ra2,dec2,sep = numpy.loadtxt('cosmos_u%d_0_8arc.txt'%(i+1),unpack=True, usecols=[0,1,-1])
 ra_2 =ra + ((-0.041 * ra + 6.059)/3600.)
 dec_2=dec + ((0.058 *dec - 0.147)/3600.)
 
 ra2_2 =ra2 + ((-0.041 * ra2 + 6.059)/3600.)
 dec2_2=dec2 + ((0.058 *dec2 - 0.147)/3600.)
 
 
 print ra2[0],dec2[0]
 print ra2_2[0],dec2_2[0]
 
 for k in range(len(ra2)):
    index =k#random.randint(1,50)
    #print 'start'
    hdr=None
    co=None

    #continue
    
    r2,d2,o,p=w.wcs_world2pix(ra2[index],dec2[index],0,h['CRVAL1'],h['CRVAL2'])
    r2_2,d2_2,o,p=w.wcs_world2pix(ra2_2[index],dec2_2[index],0,h['CRVAL1'],h['CRVAL2'])
    if r2<=h['NAXIS1'] and d2 <=h['NAXIS2'] and r2>=0 and d2>=0:
        #print i,r2,d2,'u'
  #  else:continue
        c2 = Circle((round(r2)-1, round(d2) -1), 3, edgecolor='blue', facecolor='none')
        c4 = Circle((round(r2_2), round(d2_2) -1), 3, edgecolor='green', facecolor='none')#,
    #print c2
    #sys.exit()
        ax.add_patch(c2)
        ax.add_patch(c4)
    
    if index>=len(ra):continue
    
    r,d,o,p=w.wcs_world2pix(ra[index],dec[index],0,h['CRVAL1'],h['CRVAL2'])
    
    #r,d,o,p=w.wcs_world2pix(ra[index],dec[index],0,1.501191666680E+02,2.205833333330)
        #w2 = wcs.WCS(hdr.header)
   
    
    if r<=h['NAXIS1'] and d <=h['NAXIS2'] and r>=0 and d>=0:print i,r,d, ra[index],dec[index],'detected'
   
    c = Circle((round(r)-1, round(d) -1), 4, edgecolor='red', facecolor='none')#,
    c3 = Circle((round(r), round(d)-1 ), 4, edgecolor='red', facecolor='none')#,yellow
              #transform=ax.get_transform('fk5'))              
    #ax.add_patch(c)
    ax.add_patch(c3)
    
c2 = Circle((500, 500), 5, edgecolor='yellow', facecolor='none')
#ax.add_patch(c2)
ax.imshow(data,cmap='Greys_r')
ax.set_xlabel('Right Ascension (J2000)')
ax.set_ylabel('Declination (J2000)')
#r = Rectangle((500., 500.), 10., 10., edgecolor='yellow', facecolor='none')    
#ax.add_patch(r)
gca().invert_yaxis()
show()
    
    
    
      
    
