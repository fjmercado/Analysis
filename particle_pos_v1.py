import numpy as np
import h5py
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt
import math
from operator import truediv
############################################### Position Plot Function   ##########################################################
def pos_plot(iPartType,  box_size):
    'Function plots position plots of particles for three projections'
    #... Time
    #time = isnap*0.05

    #...Tag for Saving
    tag_list = ['Gas', 'DM', '?', '',  'Stars', '' ]
    tag_type = tag_list[iPartType]

    #...Load File
    #f = h5py.File(base + '/snapshot_' + str(isnap).zfill(3) +'.hdf5')
    f = h5py.File(base + '/snapshot_555.0.hdf5')


    #...Load Stellar Data
    xyz = np.array(f['PartType'+str(iPartType)+'/Coordinates'])

    ##############################################   PLOT   ########################################################

    #...Set up figure box:
    fig = plt.figure(figsize = [20,6])


    #########################  Plot xy
    ax = fig.add_subplot(1, 3, 1)

    plt.scatter(xyz[:,0],xyz[:,1],marker='.')

    plt.axis('equal')
    plt.xlim(box_size)
    plt.ylim(box_size)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('xy')
    #ax.text(-lim+1, lim+1, str(time)+'Gyr', fontsize=15)
    
    #########################  Plot xz

    ax = fig.add_subplot(1, 3, 2)
    
    plt.scatter(xyz[:,0],xyz[:,2],marker='.')

    plt.axis('equal')
    plt.xlim(box_size)
    plt.ylim(box_size)
    plt.xlabel('x')
    plt.ylabel('z')
    plt.title('xz')
    #ax.text(-lim+1, lim+1, str(time)+'Gyr', fontsize=15)

    #########################  Plot yz
    ax = fig.add_subplot(1, 3, 3)

    plt.scatter(xyz[:,1],xyz[:,2],marker='.')

    plt.axis('equal')
    plt.xlim(box_size)
    plt.ylim(box_size)
    plt.xlabel('y')
    plt.ylabel('z')
    plt.title('yz')
    #ax.text(-lim+1, lim+1, str(time)+'Gyr', fontsize=15)

    #################################################  SAVING ########################################################

    save_fig_file = '/data8/data/mercadf1/scratch/Plots/'+tag_type+'_pos_run_ELVIS_RomeoJuliet_snap_555_0.png'
    #save_fig_file = '/data25/rouge/gonzaa11/francisco/outputs/runmed'+str(irun)+'/Plots/pos_plots/'+tag_type+'_pos_run'+str(irun)+'_snap'+str(isnap)+'.png'
    #...Report saving:
    print "Saving : "
    print str(save_fig_file)

    #...Save Figure:
    plt.savefig(save_fig_file)

    #return 'Snapshot '+str(isnap)+' complete'

############################################### SWITCHES ##########################################################

#...Choose run:
#run = 6
PartType = 0
#snaps = 41
#...Out to radius (kpc) (KEEP AS FLOAT)
lim = 100000
box = [-lim,lim]

###############################################  BOOKKEEPING  #####################################################

#...Path to data
base = '/data8/data/mercadf1/scratch/output'
#print "---------------------------------------------------------------------------"
#...Plot stuff
#for i in range(0,snaps+1):
pos_plot(PartType, box)
#   print "---------------------------------------------------------------------------"
