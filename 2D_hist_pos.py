import numpy as np
import h5py
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import math
from operator import truediv
############################################### 2D hist  Function   ##########################################################
def pos_plot(isnap, irun, iPartType,  lim):
    'Function plots position plots of particles for three projections'
    #... Time
    time = isnap*0.05
    
    #...Tag for Saving
    tag_list = ['Gas', 'DM', '?', '',  'Stars', '' ]
    tag_type = tag_list[iPartType]

    #...Load File
    f = h5py.File(base + '/snapshot_' + str(isnap).zfill(3) +'.hdf5')

    #...Load Stellar Data
    xyz = np.array(f['PartType'+str(iPartType)+'/Coordinates'])
    
    ##############################################   PLOT   ########################################################

    #...Some Parameters:
    cbar = 'magma'
    nbins = 100
    cmin,cmax = 7.0, 11.0 #0, 1.5e4
    mperp = 800

    ############# Calculate counts

    counts1, xedges, yedges, img = plt.hist2d(xyz[:,0],xyz[:,1],bins=nbins,range=[[-lim,lim],[-lim,lim]])
    counts2, xedges, yedges, img = plt.hist2d(xyz[:,0],xyz[:,2],bins=nbins,range=[[-lim,lim],[-lim,lim]])
    counts3, xedges, yedges, img = plt.hist2d(xyz[:,1],xyz[:,2],bins=nbins,range=[[-lim,lim],[-lim,lim]])

    #...Set up figure box:
    plt.style.use('dark_background')
    fig = plt.figure(figsize = [20,6])

    #########################  Plot xy
    ax = fig.add_subplot(1, 3, 1)

    #counts, xedges, yedges, img = plt.hist2d(xyz[:,0],xyz[:,1],bins=nbins)
    counts1 = np.log(counts1*mperp)
    plt.imshow(counts1, origin='lower',extent=[-lim, lim, -lim, lim], cmap=cbar, vmin=cmin, vmax=cmax)
    
    ax.axis([-lim,lim,-lim,lim])
    plt.axis('equal')
    #plt.xlim(box_size)
    #plt.ylim(box_size)
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    plt.title('xy')
    #plt.colorbar()
    ax.text(-lim, lim, str(time)+'Gyr', fontsize=20)

    #########################  Plot xz

    ax = fig.add_subplot(1, 3, 2)

    #plt.hist2d(xyz[:,0],xyz[:,2],cmap=cbar,bins=nbins,weights=mass,vmin=cmin,vmax=cmax)
    
    #counts, xedges, yedges, img = plt.hist2d(xyz[:,0],xyz[:,2],bins=nbins)
    counts2 = np.log(counts2*mperp)
    plt.imshow(counts2, origin='lower',extent=[-lim, lim, -lim, lim], cmap=cbar, vmin=cmin, vmax=cmax)
    
    ax.axis([-lim,lim,-lim,lim])
    plt.axis('equal')
    #plt.xlim(box_size)
    #plt.ylim(box_size)
    plt.xlabel('x (kpc)')
    plt.ylabel('z (kpc)')
    plt.title('xz')
    #plt.colorbar()
    #ax.text(-lim+1, lim+.5, str(time)+'Gyr', fontsize=20)

    #########################  Plot yz
    ax = fig.add_subplot(1, 3, 3)

    #plt.hist2d(xyz[:,1],xyz[:,2],cmap=cbar,bins=nbins,weights=mass,vmin=cmin,vmax=cmax)

    #counts, xedges, yedges, img = plt.hist2d(xyz[:,1],xyz[:,2],bins=nbins)
    counts3 = np.log(counts3*mperp)
    plt.imshow(counts3, origin='lower',extent=[-lim, lim, -lim, lim], cmap=cbar, vmax=cmax)

    ax.axis([-lim,lim,-lim,lim])
    plt.axis('equal')
    #plt.xlim(box_size)
    #plt.ylim(box_size)
    plt.xlabel('y (kpc)')
    plt.ylabel('z (kpc)')
    plt.title('yz')
    cb = plt.colorbar()
    cb.set_label('log(Msun/bin)',fontsize = 20)
    #ax.text(-lim+1, lim+.5, str(time)+'Gyr', fontsize=20)

    #################################################  SAVING ########################################################

    save_fig_file = '/data8/data/mercadf1/output/Pegasus/high_res/aruns/run_'+ str(irun)+'/Plots/pos_plots/'+tag_type+'_2Dhist_run'+str(irun)+'_snap'+str(isnap)+'.png'
    #save_fig_file = '/data25/rouge/gonzaa11/francisco/outputs/runmed'+str(irun)+'/Plots/pos_plots/'+tag_type+'_pos_run'+str(irun)+'_snap'+str(isnap)+'.png'
    #...Report saving:
    print "Saving : "
    print str(save_fig_file)

    #...Save Figure:
    plt.savefig(save_fig_file)

    return 'Snapshot '+str(isnap)+' complete'

############################################### SWITCHES ##########################################################

#...Choose run:
run = 2
PartType = 4
snaps = 177
#...Out to radius (kpc) (KEEP AS FLOAT)
limit = 5
box = [-limit,limit]
###############################################  BOOKKEEPING  #####################################################

#...Path to data
base = '/data8/data/mercadf1/output/Pegasus/high_res/aruns/run_'+str(run) #GreenPlanet
#base = '/data25/rouge/gonzaa11/francisco/outputs/runmed'+str(run)

print "---------------------------------------------------------------------------"
#...Plot stuff
for i in range(0,snaps+1):
    pos_plot(i, run, PartType, limit)
    print "---------------------------------------------------------------------------"
