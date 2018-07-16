import numpy as np
import h5py
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt
import math
from operator import truediv

############################################### Profile Plot Function   ##########################################################
def plot_profile(isnap, irun, itype,  my_rad, my_bins):
    
    ###############################################  BOOKKEEPING  #####################################################
    sqrt_G = 2.07e-3

    type_list = ['Gas','DM','Disk','Bulge','Stars','BH']
    type_tag = type_list[itype]

    #...Path to data
    base = '/data8/data/mercadf1/output/Pegasus/high_res/aruns/run_'+str(irun)+'a'

##############################################  Histogram   ######################################################

    #...Load File
    f = h5py.File(base + '/snapshot_' + str(isnap).zfill(3) +'.hdf5','r')

    #...Load Dark Matter Data
    xyz_DM = np.array(f['PartType1/Coordinates'])
    DM_mass = np.array(f['PartType1/Masses'])

    #...Decides what to load
    if itype == 1:
        xyz = xyz_DM
        mass = DM_mass

    else:
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
    hist1, bin_edges = np.histogram(radius, bins = np.logspace(-0.5, 1.0, my_bins), range = [0,my_rad])

    #... Check part per bin
    print ''
    print 'hist1: '
    print str(hist1)
    print ''

    #...Change to array of masses in bins
    hist = hist1*800
    
    #...For Density plot
    Vol = np.zeros(len(hist))
    for j in range(0, len(hist)):
        Vol[j] = 4.0/3.0* math.pi * (bin_edges[j +1]**3 - bin_edges[j]**3)

    #...Gets rid of zeros for density
    index = hist != 0
    hist = hist[index]
    Vol = Vol[index]
    Vol[0] = 1
    Den = map(truediv, hist, Vol)


    #...Removes added edge:
    bin_edges = bin_edges[:-1]
    #bin_edges_v = bin_edges
    bin_edges = bin_edges[index]

    #...Makes profile cummulative
    #index = hist_v != 0

    hist = np.cumsum(hist)
    v_circ = sqrt_G * (hist/bin_edges)**0.5 # In units of km/sec
    
    return bin_edges, Den, v_circ


##############################################   PLOT PROFILES  ########################################################

run = 1
bins = 21
max_rad = 10.0

for i in range(0,121):
    print '-----------------------------------------------------------------------------------------------------------'
    snap = i
    time = 0.05*snap
    print 'Snapshot: '+str(snap)
    print ''

    #...Set up figure:
    fig = plt.figure(figsize = (20,10))
    fig.suptitle('High Gas LzPlus (t = '+str(time)+'Gyr)', fontsize = 35)
    x_log = 1 #... yes (1) no (0)

    ########################################## Build Profiles

    stars_x,stars_y,stars_v = plot_profile(snap, run, 4,  max_rad, bins)
    DM_x,DM_y,DM_v = plot_profile(snap, run, 1,  max_rad, bins)

    ########################################## Plot Density Profiles

    ax1 = fig.add_subplot(1, 2, 1)

    #...Plots the bin hights vs mid-point of bins:
    plt.plot(stars_x, stars_y, linewidth=2.5, color = 'm')
    plt.plot(DM_x, DM_y, linewidth=2.5, color = 'orange')

    #...Plot details
    plt.yscale('log')
    plt.xlabel(r'$r (kpc)$', fontsize = 20)
    plt.ylabel(r'$\rho(M_{\odot}/kpc)$', fontsize = 20)
    plt.title('Density', fontsize = 25)
    plt.legend(['Stars', 'Halo'], loc='upper right')

    #...Figure out axis limis
    if x_log == 1:
        plt.xscale('log')
        ax1.set_xlim([10**-1,2*10**0])
    else:
        ax1.set_xlim([0,2])

    ax1.set_ylim([1*10**3,5*10**9])

    ########################################## Dark Matter

    ax2 = fig.add_subplot(1, 2, 2)

    #...Plots the bin hights vs mid-point of bins:
    plt.plot(stars_x, stars_v, linewidth=2.5, color = 'm')
    plt.plot(DM_x, DM_v, linewidth=2.5, color = 'orange')

    #...Plot details
    #plt.yscale('log')
    plt.xlabel(r'$r (kpc)$', fontsize = 20)
    plt.ylabel(r'$V_{circ}(km/s)$', fontsize = 20)
    plt.title('V_circ', fontsize = 25)
    plt.legend(['Stars', 'Halo'], loc='upper right')

    #...Figure out axis limis
    #if x_log == 1:
    #    plt.xscale('log')
    #    ax2.set_xlim([10**-1,2*10**0])
    #else:
    ax2.set_xlim([0,2])

    ax2.set_ylim([0,30])#([2.5*10**1,2*10**7])

    #################################################  SAVING ########################################################

    save_fig_file = '/data8/data/mercadf1/output/Pegasus/high_res/aruns/run_'+str(run)+'a/Plots/Star_DM_profile_run'+str(run)+'_snap'+str(snap)+'.png'

    #...Report saving:
    print "Saving : "
    print str(save_fig_file)

    #...Save Figure:
    fig.savefig(save_fig_file)
