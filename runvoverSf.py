import numpy as np
from v_over_sigma_AGS_f import v_over_sig

#...Choose :
irun1 = 4
irun2 = 5
#####  Choose run ########
#Path to data

snap1 = np.arange(121)
snap2 = np.arange(121)

tim1 = []
Vover1 = []
elli1 = []
theta1 = []
ratio_1 = []
for i in snap1:
    out = v_over_sig(base_path,run[0],i)
    #print out
    t=(i)*0.05
    tim1.append(t)
    ratio=out[0]
    Vover1.append(ratio)
    ell = out[1]
    elli1.append(ell)
    #rat = out[2]
    #ratio_1.append(rat)


# just in case, dealing with NaNs    
ind1 = np.where(~np.isnan(Vover1))[0]
time1 = np.array(tim1)[ind1]
VoverS1 = np.array(Vover1)[ind1]
ellip1 = np.array(elli1)[ind1]
#Ratio1 = np.array(ratio_1)[ind1]

tim2 = []
Vover2 = []
elli2 = []
theta2 = []
ratio_2 = []
for i in snap2:
    out = v_over_sig(base_path2,run[1],i)
    #print out
    t=(i)*0.05
    tim2.append(t)
    ratio=out[0]
    Vover2.append(ratio)
    ell = out[1]
    elli2.append(ell)
    #rat = out[2]
    #ratio_2.append(rat)


# just in case, dealing with NaNs    
ind2 = np.where(~np.isnan(Vover2))[0]
time2 = np.array(tim2)[ind2]
VoverS2 = np.array(Vover2)[ind2]
ellip2 = np.array(elli2)[ind2]
#Ratio2 = np.array(ratio_2)[ind2]

tim7 = []
Vover7 = []
elli7 = []
theta7 = []
ratio_7 = []
for i in snap7:
    out = v_over_sig(base_path2,run[6],i)
    t=(i)*0.05
    tim7.append(t)
    ratio=out[0]
    Vover7.append(ratio)
    ell = out[1]
    elli7.append(ell)
    #rat = out[2]
    #ratio_7.append(rat)
    #the = out[2]
    #theta7.append(the)
    
ind7 = np.where(~np.isnan(Vover7))[0]
time7 = np.array(tim7)[ind7]
VoverS7 = np.array(Vover7)[ind7]
ellip7 = np.array(elli7)[ind7]
#Ratio7 = np.array(ratio_7)[ind7]

tim8 = []
Vover8 = []
elli8 = []
theta8 = []
ratio_8 = []
for i in snap8:
    out = v_over_sig(base_path2,run[7],i)
    t=(i)*0.05
    tim8.append(t)
    ratio=out[0]
    Vover8.append(ratio)
    ell = out[1]
    elli8.append(ell)
    #rat = out[2]
    #ratio_8.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind8 = np.where(~np.isnan(Vover8))[0]
time8 = np.array(tim8)[ind8]
VoverS8 = np.array(Vover8)[ind8]
ellip8 = np.array(elli8)[ind8]
#Ratio8 = np.array(ratio_8)[ind8]

tim10 = []
Vover10 = []
elli10 = []
theta10 = []
ratio_10 = []
for i in snap10:
    out = v_over_sig(base_path2,run[9],i)
    t=(i)*0.05
    tim10.append(t)
    ratio=out[0]
    Vover10.append(ratio)
    ell = out[1]
    elli10.append(ell)
    #rat = out[2]
    #ratio_10.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind10 = np.where(~np.isnan(Vover10))[0]
time10 = np.array(tim10)[ind10]
VoverS10 = np.array(Vover10)[ind10]
ellip10 = np.array(elli10)[ind10]
#Ratio10 = np.array(ratio_10)[ind10]

tim11 = []
Vover11 = []
elli11 = []
theta11 = []
ratio_11 = []
for i in snap11:
    out = v_over_sig(base_path2,run[10],i)
    t=(i)*0.05
    tim11.append(t)
    ratio=out[0]
    Vover11.append(ratio)
    ell = out[1]
    elli11.append(ell)
    #rat = out[2]
    #ratio_11.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind11 = np.where(~np.isnan(Vover11))[0]
time11 = np.array(tim11)[ind11]
VoverS11 = np.array(Vover11)[ind11]
ellip11 = np.array(elli11)[ind11]
#Ratio11 = np.array(ratio_11)[ind11]

tim12 = []
Vover12 = []
elli12 = []
theta12 = []
ratio_12 = []
for i in snap12:
    out = v_over_sig(base_path2,run[11],i)
    t=(i)*0.05
    tim12.append(t)
    ratio=out[0]
    Vover12.append(ratio)
    ell = out[1]
    elli12.append(ell)
    #rat = out[2]
    #ratio_12.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind12 = np.where(~np.isnan(Vover12))[0]
time12 = np.array(tim12)[ind12]
VoverS12 = np.array(Vover12)[ind12]
ellip12 = np.array(elli12)[ind12]
#Ratio12 = np.array(ratio_12)[ind12]

tim13 = []
Vover13 = []
elli13 = []
theta13 = []
ratio_13 = []
for i in snap13:
    out = v_over_sig(base_path2,run[12],i)
    t=(i)*0.05
    tim13.append(t)
    ratio=out[0]
    Vover13.append(ratio)
    ell = out[1]
    elli13.append(ell)
    #rat = out[2]
    #ratio_13.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind13 = np.where(~np.isnan(Vover13))[0]
time13 = np.array(tim13)[ind13]
VoverS13 = np.array(Vover13)[ind13]
ellip13 = np.array(elli13)[ind13]
#Ratio13 = np.array(ratio_13)[ind13]


tim14 = []
Vover14 = []
elli14 = []
theta14 = []
ratio_14 = []
for i in snap14:
    out = v_over_sig(base_path2,run[13],i)
    t=(i)*0.05
    tim14.append(t)
    ratio=out[0]
    Vover14.append(ratio)
    ell = out[1]
    elli14.append(ell)
    #rat = out[2]
    #ratio_14.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind14 = np.where(~np.isnan(Vover14))[0]
time14 = np.array(tim14)[ind14]
VoverS14 = np.array(Vover14)[ind14]
ellip14 = np.array(elli14)[ind14]
#Ratio14 = np.array(ratio_14)[ind14]

tim15 = []
Vover15 = []
elli15 = []
theta15 = []
ratio_15 = []
for i in snap15:
    out = v_over_sig(base_path2,run[14],i)
    t=(i)*0.05
    tim15.append(t)
    ratio=out[0]
    Vover15.append(ratio)
    ell = out[1]
    elli15.append(ell)
    #rat = out[2]
    #ratio_15.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind15 = np.where(~np.isnan(Vover15))[0]
time15 = np.array(tim15)[ind15]
VoverS15 = np.array(Vover15)[ind15]
ellip15 = np.array(elli15)[ind15]
#Ratio15 = np.array(ratio_15)[ind15]

tim16 = []
Vover16 = []
elli16 = []
theta16 = []
ratio_16 = []
for i in snap16:
    out = v_over_sig(base_path2,run[15],i)
    t=(i)*0.05
    tim16.append(t)
    ratio=out[0]
    Vover16.append(ratio)
    ell = out[1]
    elli16.append(ell)
    #rat = out[2]
    #ratio_16.append(rat)
    #the = out[2]
    #theta8.append(the)
    
ind16 = np.where(~np.isnan(Vover16))[0]
time16 = np.array(tim16)[ind16]
VoverS16 = np.array(Vover16)[ind16]
ellip16 = np.array(elli16)[ind16]
#Ratio16 = np.array(ratio_16)[ind16]


# For plotting
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from pylab import *
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as co

mpl.rcParams['lines.linewidth']=2.5#4
mpl.rcParams['axes.linewidth']=2.0#3
mpl.rcParams['lines.markersize']=5
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.major.pad'] = 5
mpl.rcParams['xtick.minor.size'] = 7
mpl.rcParams['xtick.minor.width'] = 3
mpl.rcParams['xtick.minor.pad'] = 5
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.major.pad'] = 5
mpl.rcParams['ytick.minor.size'] = 7
mpl.rcParams['ytick.minor.width'] = 3
mpl.rcParams['ytick.minor.pad'] = 5
mpl.rcParams['font.size']=20
mpl.rcParams['xtick.labelsize']=15
mpl.rcParams['ytick.labelsize']=15

font = {'size'   : 16}


fig, (ax1) = plt.subplots(1,1, sharex=True)
fig_ratio = 2.2
subplots_adjust(hspace=0.001)
#subplots_adjust(wspace=0.001)


#### PANEL 1 #####################
ax1 = subplot(111)
ax1.set_rasterization_zorder(1)
#ax1.tick_params(labelbottom="off")

minorLocatorx = MultipleLocator(0.1)
minorLocatory = MultipleLocator(0.1)
ax1.xaxis.set_minor_locator(minorLocatorx)
ax1.yaxis.set_minor_locator(minorLocatory)
yticks = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
plt.yticks(yticks,yticks)
ax1=plt.plot(time1,VoverS1,linestyle='--',color='black',zorder=5,label=r'$King-Spheric$')
#ax1=plt.plot(time2,VoverS2,linestyle='--',color='brown',zorder=5,label=r'$Hernquist-Spheric$')
#ax1=plt.plot(time11,VoverS11,linestyle='--',color='magenta',zorder=5,label=r'$King-Spheric-Switch$')
#ax1=plt.plot(time1,Ratio1,linestyle=':',color='black',zorder=5,label=r'$King-Spheric$')
#ax1=plt.plot(time2,Ratio2,linestyle=':',color='brown',zorder=5,label=r'$Hernquist-Spheric$')
#ax1=plt.plot(time4,VoverS4,linestyle='-',color='blue',zorder=5,label='run4')
#ax1=plt.plot(time5,VoverS5,linestyle='-',color='red',zorder=5,label='run5')
#ax1=plt.plot(time6,VoverS6,linestyle='-',color='cyan',zorder=5,label='run6')
#ax1=plt.plot(time7,VoverS7,linestyle='-',color='black',zorder=5,label=r'$King-L_{z}^{+}$')
ax1=plt.plot(time12,VoverS12,linestyle='-',color='red',zorder=5,label=r'$M_{*} = 1.1\times 10^{7} M_{\odot}$')#label=r'$King-L_{z}^{+}$')
ax1=plt.plot(time13,VoverS13,linestyle='-',color='blue',zorder=5,label=r'$M_{*} = 6.9\times 10^{6} M_{\odot}$')#label=r'$King-L_{z}^{+}-high$')
ax1=plt.plot(time15,VoverS15,linestyle='-',color='magenta',zorder=5,label=r'$M_{*} = 4.2\times 10^{6} M_{\odot}$')
ax1=plt.plot(time10,VoverS10,linestyle='-',color='cyan',zorder=5,label=r'$M_{*} = 3.5\times 10^{6} M_{\odot}$')

#ax1=plt.plot(time14,VoverS14,linestyle='-',color='cyan',zorder=5,label=r'$King-L_{z}^{+}-high2$')

ax1=plt.plot(time16,VoverS16,linestyle='-',color='lime',zorder=5,label=r'$M_{*} = 1.8\times 10^{6} M_{\odot}$')
#ax1=plt.plot(time8,VoverS8,linestyle='-',color='brown',zorder=5,label=r'$Hernquist-L_{z}^{+}$')
#ax1=plt.plot(time7,Ratio7,linestyle=':',color='black',zorder=5,label=r'$King-L_{z}^{+}$')
#ax1=plt.plot(time8,Ratio8,linestyle=':',color='brown',zorder=5,label=r'$Hernquist-L_{z}^{+}$')
#ax1=plt.plot(time9,VoverS9,linestyle='-',color='grey',zorder=5,label='run9')
#ax1=plt.plot(time10,VoverS10,linestyle='-',color='lime',zorder=5,label='run10')
#ax1=plt.plot(time12,VoverS12,linestyle='-',color='green',zorder=5,label='run12')
plt.axhline(y=1.7, xmin=0.0, xmax=0.3,color='lightblue') # Aquarius
plt.axhline(y=1.99, xmin=0.0, xmax=0.3,color='lightblue') # LeoA
plt.axhline(y=1.43, xmin=0.0, xmax=0.3,color='lightblue') # Pegasus
plt.axhline(y=1.01, xmin=0.0, xmax=0.3,color='lightblue') # WLM
plt.text(0.35, 1.65, r'Aquarius ', size=14)#fontdict=font)
plt.text(0.35, 1.95, r'LeoA ', size=14)
plt.text(0.35, 1.38, r'Pegaus ', size=14)
plt.text(0.35, 0.95, r'WLM ', size=14)

ax1=plt.legend(loc= 0,prop={'size':17},title='')
#ax1=plt.title(r'$2\times V_{z}$')
plt.ylim(0.,2.2)
plt.ylabel(r'$V_{rot} / \sigma_{*} $',size=20)


plt.xlabel(' t [Gyr]',size=20)

fig.show()







##fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9,3, sharex=True)
#fig, (ax1, ax2) = plt.subplots(2,2, sharex=True)
#fig_ratio = 2.2
#subplots_adjust(hspace=0.001)
##subplots_adjust(wspace=0.001)
#
#
##### PANEL 1 #####################
#ax1 = subplot(211)
#ax1.set_rasterization_zorder(1)
#ax1.tick_params(labelbottom="off")
#
#minorLocatorx = MultipleLocator(0.5)
##minorLocatory = MultipleLocator(0.1)
#ax1.xaxis.set_minor_locator(minorLocatorx)
##ax1.yaxis.set_minor_locator(minorLocatory)
#
#ax1=plt.plot(time1,VoverS1,linestyle='--',color='black',zorder=5,label=r'$King-Spheric$')
#ax1=plt.plot(time2,VoverS2,linestyle='--',color='brown',zorder=5,label=r'$Hernquist-Spheric$')
##ax1=plt.plot(time1,Ratio1,linestyle=':',color='black',zorder=5,label=r'$King-Spheric$')
##ax1=plt.plot(time2,Ratio2,linestyle=':',color='brown',zorder=5,label=r'$Hernquist-Spheric$')
##ax1=plt.plot(time4,VoverS4,linestyle='-',color='blue',zorder=5,label='run4')
##ax1=plt.plot(time5,VoverS5,linestyle='-',color='red',zorder=5,label='run5')
##ax1=plt.plot(time6,VoverS6,linestyle='-',color='cyan',zorder=5,label='run6')
#ax1=plt.plot(time7,VoverS7,linestyle='-',color='black',zorder=5,label=r'$King-L_{z}^{+}$')
#ax1=plt.plot(time8,VoverS8,linestyle='-',color='brown',zorder=5,label=r'$Hernquist-L_{z}^{+}$')
##ax1=plt.plot(time7,Ratio7,linestyle=':',color='black',zorder=5,label=r'$King-L_{z}^{+}$')
##ax1=plt.plot(time8,Ratio8,linestyle=':',color='brown',zorder=5,label=r'$Hernquist-L_{z}^{+}$')
##ax1=plt.plot(time9,VoverS9,linestyle='-',color='grey',zorder=5,label='run9')
##ax1=plt.plot(time10,VoverS10,linestyle='-',color='lime',zorder=5,label='run10')
##ax1=plt.plot(time12,VoverS12,linestyle='-',color='green',zorder=5,label='run12')
#
#ax1=plt.legend(loc= 0,prop={'size':15},title='')
##ax1=plt.title(r'$2\times V_{z}$')
#plt.ylim(0.,1.2)
#plt.ylabel(r'$V_{rot} / \sigma_{*} $',size=20)
#
#ax2 = subplot(212)
#
#minorLocatorx = MultipleLocator(0.5)
##minorLocatory = MultipleLocator(0.1)
#ax2.xaxis.set_minor_locator(minorLocatorx)
##ax2.yaxis.set_minor_locator(minorLocatory)
#ax2=plt.plot(time1,ellip1,linestyle='--',color='black',zorder=5,label=r'$King-Spheric$')
##ax2=plt.plot(time4,ellip4,linestyle='-',color='black',zorder=5,label='run4')
##ax2=plt.plot(time5,ellip5,linestyle='-',color='red',zorder=5,label='run5')
##ax2=plt.plot(time6,ellip6,linestyle='-',color='cyan',zorder=5,label='run6')
#ax2=plt.plot(time7,ellip7,linestyle='-',color='black',zorder=5,label=r'$King-L_{z}^{+}$')
#ax2=plt.plot(time8,ellip8,linestyle='-',color='blue',zorder=5,label=r'$Hernquist-L_{z}^{+}$')
##ax2=plt.plot(time9,ellip9,linestyle='-',color='grey',zorder=5,label='run9')
##ax2=plt.plot(time10,ellip10,linestyle='-',color='lime',zorder=5,label='run10')
##ax2=plt.plot(time12,ellip12,linestyle='-',color='green',zorder=5,label='run10')
##ax1=plt.plot(time7,VoverS7,linestyle='-',color='magenta',zorder=5,label='run7')
#plt.legend(loc= 0,prop={'size':15})
#
#plt.ylabel(r'$\epsilon $',size=20)
##matplotlib.rc('ytick', labelsize=20)
##plt.xlim(7.3,11)
#plt.ylim(0.,1.)
#plt.xlabel(' t [Gyr]',size=20)
#
#fig.show()
