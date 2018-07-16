'''
Let's make  a four panel figure with the following configureation of plots:
    
    Top Left: Density profile (rho vs R) up to 5kpc [No Physics/No Feedback]
    Top RIght: Density profile (rho vs R) up to 5kpc [With Physics/With Feedback]
    Bottom Left: Vcir vs R up to 5kpc [No Physics/No Feedback]
    Bottom Right: Vcir vs R up to 5kpc [With Physics/With Feedback]
    
Make for separate components (DM, Gas, Stars)
'''
import numpy as np
import h5py
import matplotlib.pyplot as plt
import math
from operator import truediv

############################################### SWITCHES ##########################################################

#...Choose particle type (0=gas, 1=DM, 4=stars, 5=BH):
itype = 1

#...Choose run:
irun_off = 2
irun_on = 3

#...Number of bins in histogram 
my_bins = 50

#...Out to radius (kpc) (KEEP AS FLOAT)
my_rad = 5.0
    
sqrt_G = 2.07e-3 

###############################################  BOOKKEEPING  #####################################################

#... For naming file
type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
type_tag = type_list[itype]

#...Path to data
base_off = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun_off) #GreenPlanet
base_on = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun_on) #GreenPlanet

##############################################  Histograms: Physics & Feedback off   ######################################################
    
for i in range(0,151,50):
    
    #...Repot Snapshot
    print "---------------------------------------------------------------------------------------------"
    print "Snapshot: "+str(i)
    print ""
    
    #...Load File
    f_off = h5py.File(base_off + 'snapshot_' + str(i).zfill(3) +'.hdf5')
    
    #...Load Dark Matter Data
    xyz_DM = np.array(f_off['PartType1/Coordinates'])
    DM_mass = np.array(f_off['PartType1/Masses'])
    
    #...Decides what to load
    if itype == 1:
        xyz = xyz_DM
        mass = DM_mass
    
    if itype != 1:
    
        #...Load Particle Data
        xyz = np.array(f_off['PartType'+str(itype)+'/Coordinates'])
        mass = np.array(f_off['PartType'+str(itype)+'/Masses'])

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

    #...Multiply by 10**10
    mass = (mass * 10**10)
    
    #...Calculate profile
    hist, bin_edges = np.histogram(radius, bins = my_bins, weights = mass, range = [0,my_rad])

    #...Removes added edge:
    bin_edges = bin_edges[1:]

    ###############...For Density plot
    Vol = np.zeros(len(hist))
    for j in range(0, len(hist)):
        Vol[j] = 4.0/3.0* math.pi * (bin_edges[j +1]**3 - bin_edges[j]**3) 
    
    #...Gets rid of zeros for density
    index = hist != 0
    hist = hist[index]
    Vol = Vol[index] 
    Vol[0] = 1
        
    Den = map(truediv, hist, Vol)
    
    ###############...For V_cir plot
    #...Makes profile cummulative
    hist = np.cumsum(hist)
    
    v_circ = sqrt_G * (hist/bin_edges)**0.5 # In units of km/sec

#...For plotting
    if i == 0:
        bin_edges0 = bin_edges
        Den0 = Den
        v_circ0 = v_circ
        print "Snapshot 0 ready for plotting."
        
        
    if i == 50:
        bin_edges1 = bin_edges
        Den1 = Den
        v_circ1 = v_circ
        print "Snapshot 50 ready for plotting."
        
    if i == 100:
        bin_edges2 = bin_edges
        Den2 = Den
        v_circ2 = v_circ
        print "Snapshot 100 ready for plotting."
        
        
    if i == 150:
        bin_edges3 = bin_edges
        Den3 = Den
        v_circ3 = v_circ
        print "Snapshot 150 ready for plotting."
        
        
        
##############################################  Histograms: Physics & Feedback on   ######################################################
    
for i in range(0,101,50):
    
    #...Repot Snapshot
    print "---------------------------------------------------------------------------------------------"
    print "Snapshot: "+str(i)
    print ""
    
    #...Load File
    f_on = h5py.File(base_on + 'snapshot_' + str(i).zfill(3) +'.hdf5')
    
    #...Load Dark Matter Data
    xyz_DM = np.array(f_on['PartType1/Coordinates'])
    DM_mass = np.array(f_on['PartType1/Masses'])
    
    #...Decides what to load
    if itype == 1:
        xyz = xyz_DM
        mass = DM_mass
    
    if itype != 1:
    
        #...Load Particle Data
        xyz = np.array(f_on['PartType'+str(itype)+'/Coordinates'])
        mass = np.array(f_on['PartType'+str(itype)+'/Masses'])

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

    #...Multiply by 10**10
    mass = (mass * 10**10)
    
    #...Calculate profile
    hist, bin_edges = np.histogram(radius, bins = my_bins, weights = mass, range = [0,my_rad])

    #...Removes added edge:
    bin_edges = bin_edges[1:]

    ###############...For Density plot
    Vol = np.zeros(len(hist))
    for j in range(0, len(hist)):
        Vol[j] = 4.0/3.0* math.pi * (bin_edges[j +1]**3 - bin_edges[j]**3) 
    
    #...Gets rid of zeros for density
    index = hist != 0
    hist = hist[index]
    Vol = Vol[index] 
    Vol[0] = 1
        
    Den = map(truediv, hist, Vol)
    
    ###############...For V_cir plot
    #...Makes profile cummulative
    hist = np.cumsum(hist)
    
    v_circ = sqrt_G * (hist/bin_edges)**0.5 # In units of km/sec

#...For plotting
    if i == 0:
        bin_edges4 = bin_edges
        Den4 = Den
        v_circ4 = v_circ
        print "Snapshot 0 ready for plotting."
        
        
    if i == 50:
        bin_edges5 = bin_edges
        Den5 = Den
        v_circ5 = v_circ
        print "Snapshot 50 ready for plotting."
        
    if i == 100:
        bin_edges6 = bin_edges
        Den6 = Den
        v_circ6 = v_circ
        print "Snapshot 100 ready for plotting."
        
        
    if i == 150:
        bin_edges7 = bin_edges
        Den7 = Den
        v_circ7 = v_circ
        print "Snapshot 150 ready for plotting."




##############################################   PLOT PROFILES  ########################################################

#...Set up figure box:
fig = plt.figure(figsize = (14,12))


#################################  Density No Physics

ax = fig.add_subplot(2, 2, 1)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges0, Den0, linewidth=2.5, color = 'c')
plt.plot(bin_edges1, Den1, linewidth=2.5, color = 'y')
plt.plot(bin_edges2, Den2, linewidth=2.5, color = 'r')
plt.plot(bin_edges3, Den3, linewidth=2.5, color= 'k')

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
plt.xlabel(r'$radius\,[kpc]$', fontsize=15)

if itype == 4:
    plt.ylabel(r'$\rho_*\,[M_{\odot}\,kpc^{-3}]$', fontsize=15)
else:
    plt.ylabel(r'$\rho\,[M_{\odot}\,kpc^{-3}]$', fontsize=15)


#################################  Density with Physics

ax = fig.add_subplot(2, 2, 2)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges4, Den4, linewidth=2.5, color = 'c')
plt.plot(bin_edges5, Den5, linewidth=2.5, color = 'y')
plt.plot(bin_edges6, Den6, linewidth=2.5, color = 'r')

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
plt.xlabel(r'$radius\,[kpc]$', fontsize=15)

if itype == 4:
    plt.ylabel(r'$\rho_*\,[M_{\odot}\,kpc^{-3}]$', fontsize=15)
else:
    plt.ylabel(r'$\rho\,[M_{\odot}\,kpc^{-3}]$', fontsize=15)


#################################  V Circ No Physics

ax = fig.add_subplot(2, 2, 3)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges0, v_circ0, linewidth=2.5, color = 'c')
plt.plot(bin_edges1, v_circ1, linewidth=2.5, color = 'y')
plt.plot(bin_edges2, v_circ2, linewidth=2.5, color = 'r')
plt.plot(bin_edges3, v_circ3, linewidth=2.5, color= 'k')

#...Set up axis labels:
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)
plt.ylabel(r'$V_{cir}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='top right')


#################################  TOTAL Profile

ax = fig.add_subplot(2, 2, 4)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges4, v_circ4, linewidth=2.5, color = 'c')
plt.plot(bin_edges5, v_circ5, linewidth=2.5, color = 'y')
plt.plot(bin_edges6, v_circ6, linewidth=2.5, color = 'r')

#...Set up axis labels:
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)
plt.ylabel(r'$V_{cir}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='lower right')



#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun_on)+'/Plots/Pegasus_vcir.png' #Greenplanet

#save_fig_file = '/Users/Francisco/DATA/Pegasus/run_'+ str(irun)+'/plots/Pegasus_vcir_run'+str(irun)+'.png' #Local

    
#...Report saving:
print "Saving : "
print str(save_fig_file)
    
#...Save Figure:
fig.savefig(save_fig_file)





