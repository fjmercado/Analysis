import numpy as np
import h5py
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt
import math
from operator import truediv

############################################### SWITCHES ##########################################################

#...Choose run:
irun = 3

#...Number of bins in histogram 
my_bins = 50

#...Out to radius (kpc) (KEEP AS FLOAT)
my_rad = 100
    
sqrt_G = 2.07e-3 

###############################################  BOOKKEEPING  #####################################################

#...Path to data
base = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun) #GreenPlanet
#base = '/Users/Francisco/DATA/Pegasus/run_'+ str(irun) #Local

##############################################  Histogram   ######################################################
    
for i in range(0,101,50):
    
    #...Repot Snapshot
    print "---------------------------------------------------------------------------------------------"
    print "Snapshot: "+str(i)
    print ""
    
    #...Load File
    f = h5py.File(base + '/snapshot_' + str(i).zfill(3) +'.hdf5')

###############... Load Data  ###########################
    
    #...Load Dark Matter Data
    xyz_DM = np.array(f['PartType1/Coordinates'])
    DM_mass = np.array(f['PartType1/Masses'])
    
    #...Load Gas Data
    xyz_gas = np.array(f['PartType0/Coordinates'])
    gas_mass = np.array(f['PartType0/Masses'])
    
    #...Load Stellar Data
    xyz_star = np.array(f['PartType4/Coordinates'])
    star_mass = np.array(f['PartType4/Masses'])
    
########################################################
    
    #...Center of Mass Calculation (DM Center of Mass)
    DM_mass_tot = DM_mass[0]*len(DM_mass)
    mp = 0
    for j in range (0, len(DM_mass)):
        mp = DM_mass[j] * xyz_DM[j,:] + mp
    DM_cm = mp/DM_mass_tot
    print "      Center of mass (DM): " +str(DM_cm)
    
    #...Recenter:
    xyz_gas = xyz_gas - DM_cm
    xyz_star = xyz_star - DM_cm

    #...Creates an array of particle distance from center:
    radius_DM = np.sqrt(xyz_DM[:,0]**2 + xyz_DM[:,1]**2 + xyz_DM[:,2]**2)
    radius_gas = np.sqrt(xyz_gas[:,0]**2 + xyz_gas[:,1]**2 + xyz_gas[:,2]**2)
    radius_star = np.sqrt(xyz_star[:,0]**2 + xyz_star[:,1]**2 + xyz_star[:,2]**2)

    #...Get rid of hubble param    
    radius_DM = radius_DM * 0.7
    radius_gas = radius_gas * 0.7
    radius_star = radius_star * 0.7

    #...mult 10**10
    DM_mass = (DM_mass * 10 **10)
    gas_mass = (gas_mass * 10 **10)
    star_mass = (star_mass * 10 **10)

#################################  DM Profile
    #...Calculate profile
    hist_DM, bin_edges_DM = np.histogram(radius_DM, bins = my_bins, weights = DM_mass, range = [0,my_rad])
    
    #...Makes profile cummulative
    hist_DM = np.cumsum(hist_DM)
    
    #...Removes added edge:
    bin_edges = bin_edges_DM[1:] # Use this for the rest
    
    v_circ_DM = sqrt_G * (hist_DM/bin_edges)**0.5 # In units of km/sec
    
#################################  GAS Profile   
    #...Calculate profile
    hist_gas, bin_edges_gas = np.histogram(radius_gas, bins = my_bins, weights = gas_mass, range = [0,my_rad])
    
    #...Makes profile cummulative
    hist_gas = np.cumsum(hist_gas)
    
    v_circ_gas = sqrt_G * (hist_gas/bin_edges)**0.5 # In units of km/sec
    
#################################  STARS  Profile 
    #...Calculate profile
    hist_star, bin_edges_star = np.histogram(radius_star, bins = my_bins, weights = star_mass, range = [0,my_rad])
    
    #...Makes profile cummulative
    hist_star = np.cumsum(hist_star)
    
    v_circ_star = sqrt_G * (hist_star/bin_edges)**0.5 # In units of km/sec
    
#################################  TOTAL Mass profile
    
    #...Add the binned masses 
    tot_mass = hist_DM + hist_gas + hist_star
    
    v_circ_tot = sqrt_G * (tot_mass/bin_edges)**0.5 # In units of km/sec
   
    #...For plotting
    if i == 0:
        v_circ_DM0 = v_circ_DM
        v_circ_gas0 = v_circ_gas
        v_circ_star0 = v_circ_star
        v_circ_tot0 = v_circ_tot
        print "Snapshot 0 ready for plotting."
        
    if i == 50:
        v_circ_DM1 = v_circ_DM
        v_circ_gas1 = v_circ_gas
        v_circ_star1 = v_circ_star
        v_circ_tot1 = v_circ_tot
        print "Snapshot 50 ready for plotting."
        
    if i == 100:
        v_circ_DM2 = v_circ_DM
        v_circ_gas2 = v_circ_gas
        v_circ_star2 = v_circ_star
        v_circ_tot2 = v_circ_tot
        print "Snapshot 100 ready for plotting."
        
    if i == 150:
        v_circ_DM3 = v_circ_DM
        v_circ_gas3 = v_circ_gas
        v_circ_star3 = v_circ_star
        v_circ_tot3 = v_circ_tot
        print "Snapshot 150 ready for plotting."
    
        
##############################################   PLOT PROFILES  ########################################################

#...Set up figure box:
fig = plt.figure(figsize = (12,10))

#################################  DM Profile

ax = fig.add_subplot(2, 2, 1)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges, v_circ_DM0, linewidth=2.5, color = 'c')
plt.plot(bin_edges, v_circ_DM1, linewidth=2.5, color = 'y')
plt.plot(bin_edges, v_circ_DM2, linewidth=2.5, color = 'r')
#plt.plot(bin_edges, v_circ_DM3, linewidth=2.5, color= 'k')

#...Set up axis labels:
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)
plt.ylabel(r'$V_{cir,DM}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='lower right')


#################################  GAS Profile

ax = fig.add_subplot(2, 2, 2)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges, v_circ_gas0, linewidth=2.5, color = 'c')
plt.plot(bin_edges, v_circ_gas1, linewidth=2.5, color = 'y')
plt.plot(bin_edges, v_circ_gas2, linewidth=2.5, color = 'r')
#plt.plot(bin_edges, v_circ_gas3, linewidth=2.5, color= 'k')

#...Set up axis labels:
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)
plt.ylabel(r'$V_{cir,gas}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='upper right')


#################################  STAR Profile

ax = fig.add_subplot(2, 2, 3)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges, v_circ_star0, linewidth=2.5, color = 'c')
plt.plot(bin_edges, v_circ_star1, linewidth=2.5, color = 'y')
plt.plot(bin_edges, v_circ_star2, linewidth=2.5, color = 'r')
#plt.plot(bin_edges, v_circ_star3, linewidth=2.5, color= 'k')

#...Set up axis labels:
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)
plt.ylabel(r'$V_{cir,star}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='upper right')


#################################  TOTAL Profile

ax = fig.add_subplot(2, 2, 4)

#...Plots the bin hights vs mid-point of bins:
plt.plot(bin_edges, v_circ_tot0, linewidth=2.5, color = 'c')
plt.plot(bin_edges, v_circ_tot1, linewidth=2.5, color = 'y')
plt.plot(bin_edges, v_circ_tot2, linewidth=2.5, color = 'r')
#plt.plot(bin_edges, v_circ_tot3, linewidth=2.5, color= 'k')

#...Set up axis labels:
plt.xlabel(r'$radius\,[kpc]$', fontsize=18)
plt.ylabel(r'$V_{cir}[km\,sec^{-1}]$', fontsize=18)

#...Legend:
plt.legend(['IC', '0.5 Gyrs', '1 Gyrs', '1.5 Gyrs'], loc='lower right')



#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/run_'+str(irun)+'/Plots/Pegasus_vcir_r'+str(my_rad)+'.png' #Greenplanet

#save_fig_file = '/Users/Francisco/DATA/Pegasus/run_'+ str(irun)+'/plots/Pegasus_vcir_run'+str(irun)+'.png' #Local

    
#...Report saving:
print "Saving : "
print str(save_fig_file)
    
#...Save Figure:
fig.savefig(save_fig_file)
