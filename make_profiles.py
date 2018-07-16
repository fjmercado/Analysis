import numpy as np
import h5py
import matplotlib.pyplot as plt
import math
from operator import truediv

############################################### SWITCHES ##########################################################

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):
itype = 4

#...Choose run:
irun = 2

#...Number of bins in histogram 
my_bins = 40

#...Out to radius (kpc) (KEEP AS FLOAT)
my_rad = 4.

sqrt_G = 2.07e-3 

if itype == 1:
    my_rad = 56. 

###############################################  BOOKKEEPING  #####################################################

type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

#...Path to data
base = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun)

##############################################  Histogram   ######################################################
    
for i in range(0,151,50):
    
    #...Repot Snapshot
    print "---------------------------------------------------------------------------------------------"
    print "Snapshot: "+str(i)
    print ""
    
    #...Load File
    f = h5py.File(base + '/snapshot_' + str(i).zfill(3) +'.hdf5','r')
       
    #...Load Dark Matter Data
    xyz_DM = np.array(f['PartType1/Coordinates'])
    DM_mass = np.array(f['PartType1/Masses'])
    
    #...Decides what to load
    if itype == 1:
        xyz = xyz_DM
        mass = DM_mass
    
    if itype != 1:
    
        #...Load Particle Data
        xyz = np.array(f['PartType'+str(itype)+'/Coordinates'])
        mass = np.array(f['PartType'+str(itype)+'/Masses'])

    #...Center of Mass Calculation (DM Center of Mass)
    DM_mass_tot = DM_mass[0]*len(DM_mass)
    mp = 0
    for j in range (0, len(DM_mass)):
        mp = DM_mass[j] * xyz_DM[j,:] + mp
    DM_cm = mp/DM_mass_tot
    print "      Center of mass (DM): " +str(DM_cm)
    
    #...Recenter:
    xyz = xyz - DM_cm

    #...Creates an array of particle distance from center:
    radius = np.sqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2)

    #...Get rid of hubble param    
    radius = radius * 0.7
   
    #...Calculate profile
    hist, bin_edges = np.histogram(radius, bins = my_bins, weights = mass, range = [0,my_rad])
    
    #...Multiply by 10**10
    hist_use = (hist * 10**10)
    hist = (hist * 10**10)
    
    #...For Density plot
    Vol = np.zeros(len(hist))
    for j in range(0, len(hist)):
        Vol[j] = 4.0/3.0* math.pi * (bin_edges[j +1]**3 - bin_edges[j]**3) 
    
    #...Gets rid of zeros for density
    index = hist_use != 0
    hist_use = hist_use[index]
    Vol = Vol[index] 
    Vol[0] = 1
        
    Den = map(truediv, hist_use, Vol)

    #...Makes profile cummulative
    hist = np.cumsum(hist)
    
    #...Removes added edge:
    bin_edges = bin_edges[:-1]
    bin_edges_use = bin_edges
    bin_edges_use = bin_edges_use[index]
    
    v_circ = sqrt_G * (hist/bin_edges)**0.5 # In units of km/sec
    
    #...For plotting
    if i == 0:
        v_circ0 = v_circ
        bin_edges0 = bin_edges
        Den0 = Den
        bin_edges_use0 = bin_edges_use
        print "Snapshot 0 ready for plotting."
        
        
    if i == 50:
        v_circ1 = v_circ
        bin_edges1 = bin_edges
        Den1 = Den
        bin_edges_use1 = bin_edges_use
        print "Snapshot 50 ready for plotting."
        
    if i == 100:
        v_circ2 = v_circ
        bin_edges2 = bin_edges
        Den2 = Den
        bin_edges_use2 = bin_edges_use
        print "Snapshot 100 ready for plotting."
        
        
    if i == 150:
        v_circ3 = v_circ
        bin_edges3 = bin_edges
        Den3 = Den
        bin_edges_use3 = bin_edges_use
        print "Snapshot 150 ready for plotting."


##############################################   PLOT PROFILE  ########################################################
###########################... PLOT 1

#...Set up figure box:
fig = plt.figure(figsize = (10,8))
ax = fig.add_subplot(2, 1, 1)
fig.subplots_adjust(left = 0.15, right = 0.92, top = 0.89, bottom = 0.17, wspace = 0.3, hspace = 0.1)
    
#...Selects Subplot:
ax = fig.add_subplot(1)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges0, v_circ0, linewidth=2.5, color = 'c')
plt.plot(bin_edges1, v_circ1, linewidth=2.5, color = 'y')
plt.plot(bin_edges2, v_circ2, linewidth=2.5, color = 'r')
plt.plot(bin_edges3, v_circ3, linewidth=2.5, color= 'k')


if itype == 0:
    plt.ylabel(r'$V_{cir,gas}[km\,sec^{-1}]$', fontsize=18)
if itype == 1:
    plt.ylabel(r'$V_{cir,DM}[km\,sec^{-1}]$', fontsize=18)
if itype == 4:
    plt.ylabel(r'$V_{cir,star}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='lower right')

###########################... PLOT 2

ax = fig.add_subplot(2, 1, 2)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges_use0, Den0, linewidth=2.5, color = 'c')
plt.plot(bin_edges_use1, Den1, linewidth=2.5, color = 'y')
plt.plot(bin_edges_use2, Den2, linewidth=2.5, color = 'r')
plt.plot(bin_edges_use3, Den3, linewidth=2.5, color = 'k')
#plt.bar(bin_edges, hist, align='center', width=mid*2)

#...Log scale axes
plt.xscale('log')
plt.yscale('log')

#...Limit of axes:
if itype == 0:
    ax.set_xlim([10**-1,4*10**0])
    ax.set_ylim([2*10**2,4*10**8])
if itype == 1:
    ax.set_xlim([1,56])
    ax.set_ylim([10**2,2*10**7])
if itype == 4:
    ax.set_xlim([9*10**-2, 3])
    ax.set_ylim([10**2,5*10**6])


#...Set up axis labels
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)

if itype == 4:
    plt.ylabel(r'$\rho_*\,[M_{\odot}\,kpc^{-3}]$', fontsize=18)
else:
    plt.ylabel(r'$\rho\,[M_{\odot}\,kpc^{-3}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='upper right')

#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun)+'/'+type_tag+'_v_profile.png'
    
#...Report saving:
print "Saving : "
print str(save_fig_file)
    
#...Save Figure:
fig.savefig(save_fig_file)
