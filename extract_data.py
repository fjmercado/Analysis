import numpy as np
import h5py
from astropy import constants as const
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from pylab import *

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

#... Elvis (0) or LATTE (1) Andrew sims (2), and Alex sims (3)
sim = 2

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):
itype = 4

#...Snapshot
isnap = 600

#... What run
isim = 8

#... How many files per snap? 
ifile = 1

halo_id = 29843
r_half = 4.17
###############################################  BOOKKEEPING  #####################################################

#...Path to data 
if sim == 0:
    base = '/data8/data/mercadf1/scratch/output/ELVIS/'
    base1 = '/data8/data/mercadf1/scratch/Plots/ELVIS/'
    sim_list = ['nothing yet']
    
if sim == 1:
    base = '/data8/data/mercadf1/scratch/output/LATTE/'
    base1 = '/data8/data/mercadf1/scratch/Plots/LATTE/'
    sim_list = ['m12i', 'm12f', 'm12m' ]
    

if sim == 2:
    base = '/data8/data/mercadf1/scratch/output/Andrew_sims/'
    base1 = '/data8/data/mercadf1/scratch/Plots/Andrew_sims/'
    sim_list = ['h1160816', 'h194460', 'h227315', 'h286590', 'h548145', 'h548863', 'h715034', 'h738171', 'h986495']

if sim == 3:
    base = '/data8/data/mercadf1/scratch/output/Alex_sims/'
    base1 = '/data8/data/mercadf1/scratch/Plots/Alex_sims/'
    sim_list = ['halo1016', 'halo12596', 'halo20910', 'halo32503', 'halo848', 'halo948', 'halo007', 'halo11707', 'halo20192', 'halo32257', 'halo796', 'halo897']

#...File name
if sim == 0:
    fname = base + sim_list[isim] + '/snapshot_'+ str(isnap)
if sim == 1:    
    fname = base + sim_list[isim] + '/snapshot_'+ str(isnap)
if sim == 2:
    fname = base + sim_list[isim] + '/snapshot_'+ sim_list[isim] +'_Z12_bary_box_152'
if sim == 3:
    fname = base + sim_list[isim] +'_z_zero'

#...Particle Types
part_list = ['PartType0','PartType1', 'PartType2', 'PartType3','PartType4']
p_type = part_list[itype]

type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

###############################################  LOAD DATA  #####################################################
halos = np.loadtxt(base+sim_list[isim]+'/halos_0.0.ascii')
id = halos[:,0]
index = id == halo_id
x, y, z = halos[:,8], halos[:,9], halos[:,10]
x, y, z = x[index]*1000, y[index]*1000, z[index]*1000
halo_cen = np.asarray([x[0], y[0], z[0]])

#...Huble Parameter
h = get_data(fname, ifile,'Header','HubbleParam')
print 'h: '+ str(h)

#halo_cen = halo_cen

star_xyz = get_data(fname, ifile, p_type, 'Coordinates')
star_xyz = star_xyz - halo_cen

star_mass = get_data(fname, ifile, p_type, 'Masses')*(10**10)
print 'stellar data loaded...'
print str(len(star_xyz))
print str(len(star_mass))

star_z = get_data(fname, ifile, p_type, 'Metallicity')
fe_mass_frac = star_z[:,10]

r = np.sqrt(star_xyz[:,0]**2 + star_xyz[:,1]**2 + star_xyz[:,2]**2)

##############################################   Calc [Fe/H] and M_star   ########################################################

#...Everything witin r_half
index = r <= r_half 
r = r[index]
star_mass = star_mass[index]
fe_mass_frac = fe_mass_frac[index]

#...Mass (within r_half)
mass = sum(star_mass)
print 'M_star: '+str(mass)

sun_fe_frac = 0.0016

ab_ratio = np.log10(fe_mass_frac/sun_fe_frac)
print str(np.mean(ab_ratio))

##############################################   PLOT   ########################################################                                                          
#...Set up figure box:                                                                                                                                                

'''                                                                                                                                              
fig = plt.figure(figsize = (12,12))
rc('axes',linewidth=3)
plt.yticks(fontsize = 25)
plt.xticks(fontsize = 25)
plt.tick_params(which='minor',width=2,length=5)
plt.tick_params(which='major',width=2,length=10)

ax = fig.add_subplot(1,1,1)

plt.scatter(star_xyz[:,0], star_xyz[:,1], marker = '.')
#plt.scatter(halo_cen[0], halo_cen[1], marker = 'o', color = 'r')
c = plt.Circle((0,0), r_half, color = 'r', fill = False)

ax.add_patch(c)

ax.set_xlim([-10,10])                                  
ax.set_ylim([-10,10])

#################################################  SAVING ######################################################## 

save_fig_file = base1 + sim_list[isim] +'/test_plot.png'

#...Report saving:
print "Saving : "
print str(save_fig_file)

#...Save Figure:
fig.savefig(save_fig_file)
'''
