import numpy
from pylab import*
from matplotlib.ticker import AutoMinorLocator

#z1=[0.1,0.5,0.8,1.1,1.5,2.0,3.1]
#z1=[0.0,0.35,0.65,0.95,1.3,1.75,2.25,2.75,3.5,4.0]
#z1=[0.1,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.3,2.5, 2.8,3.2,3.6,4.0]
z1=[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,1.9,2.2,2.5, 2.8,3.2,3.6,4.0]
z1=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]

cosmos = numpy.loadtxt('/home/eliab/Documents/OC/Project/cosmos/cosmos_DR4.el')
def get_mag(flux):
    return -2.5*numpy.log10(flux) - 48.6
'''
surv, IDs, flx =[],[],[]
f=open('cosmos_counter.el')
lines = f.readlines()
for line in lines:
    line=line.split()
    if '#' in line:continue
    surv.append(line[3])
    IDs.append(float(line[4]))
    flx.append(float(line[8]))
surv=numpy.array(surv)
IDs=numpy.array(IDs)
flx=numpy.array(flx)    
'''

print 'done loading'
ks_lim=24.5
print 'Full sample ', len(cosmos)

gal_cosmos = cosmos[(cosmos[:,-3] - cosmos[:,-2])<0]

star_cosmos= cosmos[(cosmos[:,-3] - cosmos[:,-2])>0]


print len(gal_cosmos),'gal_cos'
print len(star_cosmos),'star_cos'
gal_cosmos2=gal_cosmos[(gal_cosmos[:,5]>=0)]
star_cosmos2=star_cosmos[(star_cosmos[:,5]>=0)]

ks_gal = get_mag(gal_cosmos[:,3])
ks_star = get_mag(star_cosmos[:,3])
hist(ks_star[(10<ks_star)&(ks_star<30)],50,histtype='step',color='blue', label='Stars')
hist(ks_gal[(10<ks_gal)&(ks_gal<30)],50,histtype='step',color='black',linestyle='dotted', label='galaxies')
xlabel('ks',fontsize=18)
#yscale('LOG')
legend()

fig,ax = subplots()
hist(star_cosmos[:,5][star_cosmos[:,5]>=0],20,histtype='step',color='blue', label='Stars')
hist(gal_cosmos[:,5][  gal_cosmos[:,5]>=0],20,histtype='step',color='black',linestyle='--', label='galaxies')
xlabel('Redshift',fontsize=18)
yscale('LOG')
legend()

ks_f_cosmos = gal_cosmos
ks_cosmos = ks_f_cosmos[(0 < get_mag(ks_f_cosmos[:,3]))& (get_mag(ks_f_cosmos[:,3])<=ks_lim)]#23.7#22.4
print len(ks_cosmos),'ks Deep'
#ks_cosmos=ks_cosmos[ks_cosmos[:,12]>0]
#s_cosmos=ks_cosmos[ks_cosmos[:,11]>0]
#print len(ks_cosmos), 'removed null z'
#ks_cosmos=ks_cosmos[ks_cosmos[:,22]>0]
#print len(ks_cosmos), 'removed null mass'

#ks= 
print   ks_cosmos[:,3]
print get_mag(ks_cosmos[:,3])
Mass_lim = ks_cosmos[:,-1] - 0.4*(ks_lim - get_mag(ks_cosmos[:,3]))

M_comp=[]
Z_med=[]
z=z1#numpy.arange(0,4.1,0.1)

f1=open('cos_data/data_s1.el','w')
f2=open('cos_data/data_s2.el','w')
f3=open('cos_data/data_s3.el','w')
f4=open('cos_data/data_s4.el','w')
f5=open('cos_data/data_s5.el','w')
f6=open('cos_data/data_s6.el','w')
f7=open('cos_data/data_s7.el','w')
f8=open('cos_data/data_s8.el','w')
f9=open('cos_data/data_s9.el','w')
f10=open('cos_data/data_s10.el','w')
f11=open('cos_data/data_s11.el','w')
f12=open('cos_data/data_s12.el','w')
f13=open('cos_data/data_s13.el','w')
f14=open('cos_data/data_s14.el','w')
f15=open('cos_data/data_s15.el','w')
f16=open('cos_data/data_s16.el','w')
f17=open('cos_data/data_s17.el','w')

f0=open('cosmos_comp.el','w')

fnames=[f1,f2,f3,f4,f5,f6,f7,f8, f9, f10, f11, f12, f13, f14, f15, f16, f17]
for i in range(len(fnames)): 
 #fl ='f%s'%i
 fnames[i].write('# %13s %15s %6s %7s \n'%('ra', 'dec', 'photo_z','ID'))

print z1
print r'\begin{tabular}{ccccc} \hline'
print r'Redshift bin & $N_\rm{Tot}$  & $\log(M_{\rm{lim}}/M_\odot)$ & $N$ & $N_\rm{VLA}$ \\ \hline'
for i in range(len(z)-1):
    #z_cosmos = ks_cosmos[(ks_cosmos[:,11] <z[i+1])&(ks_cosmos[:,11] >z[i])]
    z_cosmos = ks_cosmos[(ks_cosmos[:,5] <z[i+1])&(ks_cosmos[:,5] >z[i])]
    #Mass_lim_z = Mass_lim[(ks_cosmos[:,11] <z[i+1])&(ks_cosmos[:,11] >z[i])]
    Mass_lim_z = Mass_lim[(ks_cosmos[:,5] <z[i+1])&(ks_cosmos[:,5] >z[i])]
    m_comp = numpy.percentile(Mass_lim_z[Mass_lim_z>0],90)
    M_comp.append(m_comp)
    z_cosmos_comp=z_cosmos[z_cosmos[:,-1]>=m_comp]
    if z[i] in z1:
     #print z[i]' < z <',z[i+1], 'completeness limit ',m_comp
     print r'$',z[i],'< z <',z[i+1],r'$&$',len(z_cosmos),r'$&', round(m_comp,1),r'&',len(z_cosmos[z_cosmos[:,-1]>m_comp]), r'&'
     #print 'N_sample full',len(z_cosmos),'N_sample complete', len(z_cosmos[z_cosmos[:,22]>m_comp]) 
    else:print z[i] ,'not in z1'
    '''
    print z[i],'< z <',z[i+1]
    print 'N_sample full',len(z_cosmos)
    print 'completeness limit ',m_comp
    print 'N_sample complete', len(z_cosmos[z_cosmos[:,22]>m_comp])
    '''
    z_med=(z[i]+z[i+1])/2.
    Z_med.append(z_med)
    c=0
    for cos in z_cosmos_comp:
    #for cos in z_cosmos: #ignore mass lim
        #if cos[-1] in IDs:
         # if surv[IDs==cos[-1]] =='COSMOS2015':
          #  c+=1

        fnames[i].write('%15.10f %15.10f %6.2f %7d \n'%(cos[1],cos[2],cos[5], cos[0]))
        f0.write('%15.10f %15.10f %6.2f %7d \n'%(cos[1],cos[2],cos[5], cos[0]))
    #print c ,r' \\'
    
#print len(Mass_lim_lowz[Mass_lim_lowz<m_lim])

fig,ax = subplots()
#plot(ks_cosmos[:,5],Mass_lim,'b.',label='Mass limit', alpha=0.01)
plot(ks_cosmos[:,5],ks_cosmos[:,-1],'k.',label='Stellar mass', alpha=0.02)
plot(Z_med,M_comp, '-ro', label='Mass comp limit',fillstyle='none',markersize=10)
ylim(7.48,11.99)
xlim(0,4.9)
ylabel(r'$\log(M/M_\odot)$',fontsize=27)
xlabel(r'$z_{\rm{photo}}$',fontsize=34)
tick_params(axis='both',which = 'major', labelsize=20,width =3)
tick_params(axis='both',which = 'minor', labelsize=12, width=1)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
legend(loc='upper left')
show()

f0.close()

for i in range(len(fnames)): 
 fnames[i].close()
#sf_cosmos= ks_cosmos[ks_cosmos[:,-1]==1]
#print 'sf',len(sf_cosmos)
#ks_cosmos_lowz = sf_cosmos[(sf_cosmos[:,-2] <0.4)&(sf_cosmos[:,-2] >0.1)]
#ks_cosmos_lowz = ks_cosmos[(ks_cosmos[:,12] <0.5)&(ks_cosmos[:,12] >0.1)]
#Mass_lim_lowz = Mass_lim[(ks_cosmos[:,12] <0.5)&(ks_cosmos[:,12] >0.1)]
#print len(ks_cosmos_lowz), 'lowz'





'''
print len(ks_cosmos_lowz),'lowz'
U_lz=ks_cosmos_lowz[:,9]
V_lz=ks_cosmos_lowz[:,7]
J_lz=ks_cosmos_lowz[:,6]


gal_cosmos_lowz = gal_cosmos[(gal_cosmos[:,-2] <0.4)&(gal_cosmos[:,-2] >0.1)]

Nuv =gal_cosmos_lowz[:,-4]
print Nuv,'Nuv'

r=gal_cosmos_lowz[:,8]
print r
J = gal_cosmos_lowz[:,6]

x = V_lz - J_lz
uvj_1 = ks_cosmos_lowz[x>=1]
uvj_2 = ks_cosmos_lowz[x<1]
y_1 =0.68*(x[x>=1])+.92
y_2 =0.*(x[x<1])+ 0.68+.92
sf_uvj_2 = uvj_2[(U_lz[x<1] - V_lz[x<1])<y_2]
sf_uvj_1 = uvj_1[(U_lz[x>=1] - V_lz[x>=1])<y_1]

f1=open('data_s1.el','w')
f2=open('data_s2.el','w')
f3=open('data_s3.el','w')
f4=open('data_s4.el','w')
f5=open('data_s5.el','w')
f6=open('data_s6.el','w')
f7=open('data_s7.el','w')
f8=open('data_s8.el','w')
f9=open('data_s9.el','w')
f10=open('data_s10.el','w')
f11=open('data_s11.el','w')

fnames=[f1,f2,f3,f4,f5,f6,f7, f7,f8, f9, f10, f11]
for i in range(11): 
 #fl ='f%s'%i
 fnames[i].write('# %15s %15s %6s %6s %6s   %6s %6s %6s %6s %6s %6s %6s %6s %6s  \n'%('ra', 'dec', 'UltraV', 'peter', 'Ks','Ks_f','J','V','r','U', 'Nuv','type','z_p','sf'))
print 'sf', len(sf_cosmos)
for i in range(len(sf_cosmos)-1):
    z = sf_cosmos[i][-2]
    if z>0.1 and z<=0.4:num=0
    elif z>0.4 and z<=0.6:num=1
    elif z>0.6 and z<=0.8:num=2
    elif z>0.8 and z<=1.0:num=3
    elif z>1.0 and z<=1.3:num=4
    elif z>1.3 and z<=1.6:num=5
    elif z>1.6 and z<=2.0:num=6
    elif z>2.0 and z<=2.5:num=7
    elif z>2.5 and z<=3.3:num=8
    elif z>3.3 and z<=4.6:num=9
    elif z>4.6 and z<=5.7:num=10
    #if i ==5262:
    #   print i, sf_cosmos[i][0],sf_cosmos[i][1], len(sf_cosmos[i]),len(sf_cosmos[i-1])
    #   continue
    
    
    fnames[num].write('%15.10f %15.10f  %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f %6.2f %6.2f   %6.2f %2d %6.2f %1d \n'%(sf_cosmos[i][0],sf_cosmos[i][1],sf_cosmos[i][2],sf_cosmos[i][3],sf_cosmos[i][4],sf_cosmos[i][5] ,sf_cosmos[i][6],sf_cosmos[i][7],sf_cosmos[i][8],sf_cosmos[i][9],sf_cosmos[i][10],sf_cosmos[i][11],sf_cosmos[i][12],sf_cosmos[i][13]))

for i in range(11): 
 fnames[i].close()

f=open('uvj_lowz_2.el','w')
f.write('# %15s %15s %6s %6s %6s   %6s %6s %6s %6s %6s %6s %6s %6s %6s  \n'%('ra', 'dec', 'UltraV', 'peter', 'Ks','Ks_f','J','V','r','U', 'Nuv','type','z_p','sf'))

for i in range(len(sf_uvj_1)):

    f.write('%15.10f %15.10f  %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f %6.2f %6.2f   %6.2f %2d %6.2f %1d \n'%(sf_uvj_1[i][0],sf_uvj_1[i][1],sf_uvj_1[i][2],sf_uvj_1[i][3],sf_uvj_1[i][4],sf_uvj_1[i][5] ,sf_uvj_1[i][6],sf_uvj_1[i][7],sf_uvj_1[i][8],sf_uvj_1[i][9],sf_uvj_1[i][10],sf_uvj_1[i][11],sf_uvj_1[i][12],sf_uvj_2[i][13]))

    
for i in range(len(sf_uvj_2)):

    f.write('%15.10f %15.10f  %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f %6.2f %6.2f   %6.2f %2d %6.2f %1d \n'%(sf_uvj_2[i][0],sf_uvj_2[i][1],sf_uvj_2[i][2],sf_uvj_2[i][3],sf_uvj_2[i][4],sf_uvj_2[i][5] ,sf_uvj_2[i][6],sf_uvj_2[i][7],sf_uvj_2[i][8],sf_uvj_2[i][9],sf_uvj_2[i][10],sf_uvj_2[i][11],sf_uvj_2[i][12], sf_uvj_2[i][13]))
f.close()

fig,ax = plt.subplots()
plot((r-J),(Nuv-r),'.')
xlabel(' r - J',fontsize=25)
ylabel(' Nuv - r',fontsize=25)
tick_params(axis='both',which = 'major', labelsize=20,width =3)
tick_params(axis='both',which = 'minor', labelsize=10, width=2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
show()



#U_lz=ks_cosmos_lowz[:,11]
#V_lz=ks_cosmos_lowz[:,7]
#V_lz=ks_cosmos_lowz[:,9]
#J_lz=ks_cosmos_lowz[:,7]
#J_hz=ks_cosmos_highz[:,7]
#V_hz=ks_cosmos_highz[:,9]
#U_hz=ks_cosmos_highz[:,11]
#c_l = .69
#c_h = .59

fig,ax = plt.subplots()



plot(x[x<1],0.*(x[x<1])+ 0.68+.92,'b',linewidth=2)
plot(x[x>=1],0.68*(x[x>=1])+.92,linewidth=2)
plot((V_lz - J_lz),(U_lz - V_lz),'.k',alpha=0.3)
xlabel(' V - J',fontsize=25)
ylabel(' U - V',fontsize=25)
xlim(0,3)
ylim(0,3)
tick_params(axis='both',which = 'major', labelsize=20,width =3)
tick_params(axis='both',which = 'minor', labelsize=10, width=2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
show()
'''

