import numpy as np
import h5py
from astropy import constants as const
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from pylab import *
import math
from operator import truediv

############################################### FUNCTIONS  ##########################################################                                                                   

'''                                                                                                                                                                                      
                                                                                                                                                                                         
# filename: name of the snaps file(including the path)                                                                                                                                  
                                                                                                                                                                                         
#num_of_file: the number of the blocks                                                                                                                                                   
                                                                                                                                                                                         
#key1: 'Header','PartType0','PartType1','PartType2',',etc                                                                                                                                
                                                                                                                                                                                         
#key2: if key1=='Header' then key2 = 'Time','Redshift','HubbleParam',etc;                                                                                                                
                                                                                                                                                                                         
#      else   key2 = 'Coordinates','Masses','ParticleIDs',etc                                                                                                                            
                                                                                                                                                                                         
'''

def get_data(filename,num_of_file,key1,key2):
    if num_of_file == 1:
        f = h5py.File(filename+'.hdf5', 'r')
        if key1 == 'Header':
            return f[key1].attrs[key2]
        else:
            return f[key1][key2][:]
    else:
        for i in range(0,num_of_file):
            f = h5py.File(filename+'.'+str(i)+'.hdf5', 'r')
            if key1 == 'Header':
                return f[key1].attrs[key2]
            else:
                if ( len(f[key1][key2][:].shape)==1 ):
                    if i==0:
                        result = f[key1][key2][:]
                    else:
                        result = np.hstack( (result,f[key1][key2][:]) )
                else:
                    if i==0:
                        result = f[key1][key2][:]
                    else:
                        result = np.vstack( (result,f[key1][key2][:]) )
        return result
############################################### SWITCHES ##########################################################                                                                      
#... Elvis (0) or LATTE (1)                                                                                                                                                              
sim = 1

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):                                                                                                                                   
itype = 0

#...Snapshot                                                                                                                                                                             
isnap = 600

#... What run?                                                                                                                                                                           
irun = 'm12f'

#... How many files per snap?
ifile = 4

#...Number of bins in histogram                                                                                                                                                          
my_bins = 150

#...Out to radius (kpc) (KEEP AS FLOAT)                                                                                                                                                  
my_rad = 200.0

###############################################  BOOKKEEPING  #####################################################                                                                      

part_list = ['PartType0','PartType1','PartType2','PartType4']
p_type = part_list[itype]

type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

#...Path to data                                                                                                                                                                         

if sim == 0:
    base = '/data8/data/mercadf1/scratch/output/ELVIS'
    base1 = '/data8/data/mercadf1/scratch/Plots/ELVIS/'
if sim == 1:
    base = '/data8/data/mercadf1/scratch/output/LATTE'
    base1 = '/data8/data/mercadf1/scratch/Plots/LATTE/'


fname = base + '/' + irun + '/snapshot_'+ str(isnap)
fname1 = base + '/' + irun + '/halos_' + str(isnap)
###############################################  LOAD DATA  #####################################################

h = get_data(fname, ifile,'Header','HubbleParam')
print 'h: '+ str(h)
gas_xyz = get_data(fname, ifile, p_type, 'Coordinates')/h
gas_mass = get_data(fname, ifile, p_type, 'Masses')*(10**10)/h
print 'Gas data loaded...'
f = h5py.File(fname1 + '.hdf5', 'r')
halo_pos = np.array(f['position'])#/h#/(10**3)                                                                                                                                           
host_cent = halo_pos[0]
print 'DM data loaded...'
print "      Host halo center (DM): " +str(host_cent)

#...For Temperature calculation                                                                                                                                                          
gas_inten = get_data(fname, 4, p_type, 'InternalEnergy')
gas_z = get_data(fname, 4, p_type, 'Metallicity')
e_ab = get_data(fname, 4, p_type, 'ElectronAbundance')

############################################### CALCULATE GAS TEMP  #####################################################                                                                
#...Helium Mass fraction                                                                                                                                                                 
he_mass_frac = gas_z[:,1]
y_he = he_mass_frac / (4*(1 - he_mass_frac))

#...Mean molecular weight                                                                                                                                                                
Molecular_weights = (1.0 + 4.0*y_he)/(1.0+y_he+e_ab)*const.m_p.value
temp = (1000)**(2)*(5.0/3.0-1.0)*Molecular_weights*gas_inten / const.k_B.value #This should give temp in K                                                                               

print "Some numbers"
print "--------------------------------------------------"
print 'mean Temp all gas: ' + str(np.mean(temp))

############################################### TEMPERATURE CUT-OFF  #####################################################                                                               

index = temp >= 1e6

temp_high = temp[index]
print "Hot gas mean temp: "+ str(np.mean(temp_high))
xyz = gas_xyz[index]
mass = gas_mass[index]
print "Hot gas total mass: "+ str(np.sum(mass))

############################################### MAKE PROFILES  #####################################################                                                                     

########### DENSITY PROFILE                                                                                                                                                              

#...Recenter:                                                                                                                                                                            
xyz = xyz - host_cent

#...Creates an array of particle distance from center:                                                                                                                                   
radius = np.sqrt(xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2])

#...Calculate profile                                                                                                                                                                    
hist, bin_edges = np.histogram(radius, bins = my_bins, weights = mass, range = [0,my_rad])


#...For Density plot                                                                                                                                                                     
Vol = np.zeros(len(hist))
for j in range(0, len(hist)):
    Vol[j] = 4.0/3.0* math.pi * (bin_edges[j +1]**3 - bin_edges[j]**3)

bin_edges = bin_edges[:-1]

#...Gets rid of zeros for density                                                                                                                                                        
#index1 = hist != 0                                                                                                                                                                      
#hist = hist[index1]                                                                                                                                                                     
#Vol = Vol[index1]                                                                                                                                                                       
#bin_edges = bin_edges[index]                                                                                                                                                            
#Vol[0] = 1                                                                                                                                                                              


Den = map(truediv, hist, Vol)

bin_edges = np.asarray(bin_edges)
Den = np.asarray(Den)

print '----------------------------------'
print 'bin_edges/Den Lengths: '+ str(len(bin_edges)) + ' and ' + str(len(Den))

#...Save Profile to .csv file
prof = np.asarray([bin_edges, Den])
save_txt = base1+str(irun)+'/hot_gas_prof_snap'+str(isnap)+'.csv'
np.savetxt(save_txt, prof, delimiter=",")
print '.csv saved'
##############################################   PLOT PROFILE  ########################################################                                                                  

#...Set up figure box:                                                                                                                                                                   
fig = plt.figure(figsize = (12,12))
rc('axes',linewidth=3)
plt.yticks(fontsize = 25)
plt.xticks(fontsize = 25)
plt.tick_params(which='minor',width=2,length=5)
plt.tick_params(which='major',width=2,length=10)

ax = fig.add_subplot(1,1,1)

print bin_edges
print '--------------------------'
print Den

plt.plot(bin_edges, Den, linewidth=2.5, color = 'purple')


#...Log scale axes                                                                                                                                                                      
plt.xscale('log')
plt.yscale('log')

#ax.set_xlim([10**0,2*10**2])                                                                                                                                                            
#ax.set_ylim([1*10**-2,2*10**4])                                                                                                                                                         

plt.title(str(irun) + ' ' + str(isnap), fontsize = 30)
plt.xlabel(r'$Radius\, [kpc]$', fontsize = 30)
plt.ylabel(r'$\rho_{hot\, gas}\, [M_{\odot}\, kpc^{-3}]$', fontsize = 30)

#plt.scatter(xyz[:,0], xyz[:,1], marker = '.', c = 'r')                                                                                                                                  
#plt.scatter(host_cent[0], host_cent[1], marker = 's', c = 'y')                                                                                                                          

#################################################  SAVING ########################################################                                                                       
save_fig_file = base1+str(irun)+'/hot_gas_prof_snap'+str(isnap)+'.png'

#...Report saving:                                                                                                                                                                       
print "Saving : "
print str(save_fig_file)

#...Save Figure:                                                                                                                                                                         
fig.savefig(save_fig_file)
