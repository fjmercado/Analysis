"""
Location:     /n/home05/moreno/Work/ForStudents/Maps/codes/
Created:      M/05/03/16 (Caltech)
Last update:  M/05/03/16 (Caltech)
Author:       Jorge Moreno
Status:       Works for Moreno Runs only (GG_5 with acce)
Requirements: Access to /Runs/FIRE_Mergers/
Pipeline:     plot_map.py
Reference:    ../FaceOnProfiles/codes/plot_map.py
"""

#######################################################################

#...Import key external routines:

import numpy as np
import glob 
import sys
import h5py
import util
import simread.readsnapHDF5 as ws
import os
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm     as mplcm
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import binned_statistic

##################################################################
#                                                                #
# This creats maps of FIRE galaxy mergers.                       #
#                                                                #
##################################################################
#                                                                #
# CLIENTS: ../make_map.py.                                       #
#                                                                #
##################################################################
#                                                                #
# MEMBER FUNCTIONS:                                              #
#                                                                #
#     __init__()                                                 #
#     _setup_paths()                                             #
#     _load_bh_data()                                            #
#     _load_sim_data(comp)                                       #
#     plot_map()                                                 #
#     _setup_fig_map()                                           #
#                                                                #
##################################################################
#                                                                #
# PRODUCTS:                                                      #
#                                                                #
#    ../plots/ or ../frames/                                     #
#                                                                #
##################################################################

#...Set up font size:
fontsize = 14

#######################################################################

class Map:
    
    # __doc__:
    """This creates maps of FIRE galaxy mergers."""
    
    #...Class initialization:
    def __init__(self, iorbit = 0, isnap = 0, igal = 0, platform = 'odyssey'):

        #---Basics:
        print ""

        #...Collect basic values:
        self.platform = platform
        self.iorbit   = iorbit
        self.isnap    = isnap
        self.igal     = igal
      
        #...Declare all-inclusive data:
        self.all_data = {}
 
        #...Set up paths:
        self._setup_paths()
        
        #...Load and shuffle BH data:
        self._load_bh_data()
    
        #...Load gas data:
        self._load_sim_data(comp='gas')
        #print self.all_data['gas_pos']
    
        #...Load star data:
        #self._load_sim_data(comp='stars')
    
    ################## INITIALIZATION PRIVATE FUNCTIONS ###################
    
    #---Private Class Function: Set up paths for simulated data.
    def _setup_paths(self):
        
        #...Set up merger tag:
        self.merger_tag = 'G0'
        if self.iorbit > 0: self.merger_tag = 'GG_'+str(self.iorbit)
    
        #---On Odyssey:
        if self.platform == 'odyssey':
        
            #...Set up path to simulations:
            self.sim_path  = '/n/home05/moreno/Runs/FIRE_Mergers/'

            #...Set up path to orbit:
            self.orbit_path = 'Iso/Control_G0'
            if self.iorbit > 0: self.orbit_path = 'Merg/Merger_GG_'+str(self.iorbit)
            self.orbit_path = self.sim_path+self.orbit_path

            #...Set up path to snapshot
            self.snap_path = self.orbit_path+'/output/snapshot_'+str(self.isnap).zfill(3)+'.hdf5'
            print self.snap_path

        #---On Local Computer:
        elif self.platform == 'local':
            
            #...Set up path to simulations:
            self.user_path = '/Users/student/'
            self.snap_path = self.user_path+'Work/Mergers/Maps/hdf5_files/'+str(self.merger_tag)+'/'+str(self.merger_tag)+'_snapshot_'+str(self.isnap).zfill(3)+'.hdf5'
        print self.snap_path

    ######################################################

    #---Private Class Function: Load and shuffle BH data.
    def _load_bh_data(self):
        
        #---On Odyssey:
        if self.platform == 'odyssey':

            #...Load BH data:
            bh_id   = np.array(ws.read_block(self.snap_path,'ID  ', parttype=5))
            bh_xyz  = np.array(ws.read_block(self.snap_path,'POS ', parttype=5))
            bh_vxyz = np.array(ws.read_block(self.snap_path,'VEL ', parttype=5))
        
        #---On Local Computer:
        elif self.platform == 'local':
            
            #...Retrieve hdf5 data:
            data_hdf5 = h5py.File(str(self.snap_path),'r')
        
            #...Load BH data:
            bh_id   = np.array(data_hdf5['bh_id'])
            bh_xyz  = np.array(data_hdf5['bh_pos'])
            bh_vxyz = np.array(data_hdf5['bh_vel'])
            print "pos",bh_xyz

            print "vel",bh_vxyz

        #...Shuffle BHs:
        ord_index_bh_id = np.argsort(bh_id)

        #...Collect BH positions:
        self.all_data['bh_id']  = bh_id[ord_index_bh_id]
        self.all_data['bh_pos'] = bh_xyz[ord_index_bh_id,:]
        self.all_data['bh_vel'] = bh_vxyz[ord_index_bh_id,:]
    
        #...Compute BH separation:
        if self.iorbit != 0:
            
            dx_bh = bh_xyz[1,0]-bh_xyz[0,0]
            dy_bh = bh_xyz[1,1]-bh_xyz[0,1]
            dz_bh = bh_xyz[1,2]-bh_xyz[0,2]
            self.r_bh = np.sqrt(dx_bh*dx_bh+dy_bh*dy_bh+dz_bh*dz_bh)

    ######################################################

    #---Private Class Function: Load gas data.
    def _load_sim_data(self, comp='gas'):

        #...Set up particle type (gas or stars):
        parttype = 0
        if comp == 'stars': parttype = 4
        
        #---On Odyssey:
        if self.platform == 'odyssey':
        
            #...Load basic data:
            self.all_data[comp+'_id']   = np.array(ws.read_block(self.snap_path,'ID  ', parttype=parttype))
            self.all_data[comp+'_pos']  = np.array(ws.read_block(self.snap_path,'POS ', parttype=parttype))
            self.all_data[comp+'_vel']  = np.array(ws.read_block(self.snap_path,'VEL ', parttype=parttype))
            self.all_data[comp+'_mass'] = 1.e10*np.array(ws.read_block(self.snap_path,'MASS', parttype=parttype))
            self.all_data[comp+'_sfr']  = 1.e10*np.array(ws.read_block(self.snap_path,'SFR ', parttype=parttype))
        
            #...Load acceleration data (gas only):
            if comp == 'gas':
        
                self.all_data[comp+'_acc']       = np.array(ws.read_block(self.snap_path,'ACCE', parttype=parttype))
                self.all_data[comp+'_hydro_acc'] = np.array(ws.read_block(self.snap_path,'HACC', parttype=parttype))
                self.all_data[comp+'_grav_acc']  = self.all_data[comp+'_acc'] - self.all_data[comp+'_hydro_acc']

            #...Re-center positions:
            self.all_data[comp+'_pos']  = self.all_data[comp+'_pos'] - self.all_data['bh_pos'][self.igal,:]
            self.all_data[comp+'_vel']  = self.all_data[comp+'_vel'] - self.all_data['bh_vel'][self.igal,:]
            
            #if comp == 'gas': print "queried hdf5: ", np.asarray(self.all_data[comp+'_pos'][:,0])

        #---On Local Computer:
        elif self.platform == 'local':
        
            #...Retrieve hdf5 data:
            data_hdf5 = h5py.File(str(self.snap_path),'r')
            
            #...Load basic data:
            self.all_data[comp+'_id']   = np.array(data_hdf5[comp+'_id'])
            self.all_data[comp+'_pos']  = np.array(data_hdf5[comp+'_pos']) #JM: For your profiles
            self.all_data[comp+'_vel']  = np.array(data_hdf5[comp+'_vel']) #JM: For your profiles
            self.all_data[comp+'_mass'] = np.array(data_hdf5[comp+'_mass']) #JM: For your profiles
            self.all_data[comp+'_sfr']  = np.array(data_hdf5[comp+'_sfr']) #JM: For your profiles
            #if comp == 'gas': print "uploaded hdf5: ", np.asarray(data_hdf5[comp+'_pos'][:,0])
            #if comp == 'gas': print "copied hdf5: ", np.asarray(self.all_data[comp+'_pos'][:,0])

            #...Load acceleration data (gas only):
            if comp == 'gas':
        
                self.all_data[comp+'_acc']       = np.array(data_hdf5[comp+'_acc'])
                self.all_data[comp+'_hydro_acc'] = np.array(data_hdf5[comp+'_hydro_acc'])
                self.all_data[comp+'_grav_acc']  = self.all_data[comp+'_acc'] - self.all_data[comp+'_hydro_acc']


    ########################## PUBLIC FUNCTIONS ############################

    #---Class Function: Plots acceleration versus radius.
    def plot_map(self, comp='gas', save_tag = '', title = '', label = None, eps_tag=False, save_frame=True):
    
        #---Basic Initialization:
        
        #...Set up save tag:
        plot_tag = comp+'_map'
        save_tag = plot_tag

        #---Plot contour data:
        
        #...Set up figure:
        fig, ax = self._setup_fig_map(save_frame=save_frame)
        
        #...Set up arrays for plotting
        plot_x = self.all_data[comp+'_pos'][:,0]
        plot_y = self.all_data[comp+'_pos'][:,1]
        plot_z = self.all_data[comp+'_vel'][:,2]
#        plot_z = self.all_data[comp+'_sfr']
#        plot_z = self.all_data[comp+'_mass']
#        plot_z = self.all_data[comp+'_vel'][:,2]
        #print "plot_x:", plot_x
        
        #...Create 2D histogram:
        my_norm = None #LogNorm()
#        my_norm = colors.Normalize(vmin=-1000.,vmax=1000.)
        my_cmap = 'bwr' #'spectral' #'RdBu'#'CMRmap_r'#'ocean_r'#'cubehelix_r'
        my_cmap = copy.copy(matplotlib.cm.get_cmap(my_cmap))
        my_cmap.set_bad((0,0,0))
        my_cmap.set_bad('w',1.)
        #counts, ybins, xbins, image = ax.hist2d(plot_y, plot_x, weights=None, bins=150, range=[[-10,10],[-10,10]], cmap=my_cmap, norm=my_norm)
        counts, ybins, xbins, image = ax.hist2d(plot_y, plot_x, weights=plot_z, bins=150, range=[[-10,10],[-10,10]], cmap=my_cmap, norm=my_norm)
        ax.imshow(counts, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], interpolation='bicubic', origin='lower', norm=my_norm, cmap=my_cmap, vmin = -200., vmax = 200)
        counts = np.rot90(counts)
        counts = np.flipud(counts)
        masked_counts = np.ma.masked_where(counts==0,counts)
        ax.imshow(masked_counts, extent = [xbins.min(), xbins.max(), ybins.min(), ybins.max()], interpolation ='bicubic', origin='lower', norm=my_norm, cmap=my_cmap)


        #ColorMap Normalization For RdBu
        
        #norm(plot_z)
        #pcm = ax.pcolormesh(x, y, Z, vmin=-1., vmax=1., cmap='RdBu')
        
        #---Save contour plot:

        #...Set up basics:
        title = save_tag+'_'+self.merger_tag
        plt.title(title)

        #...Save single figure to pdf:
        if save_frame == False:

            #...If saved plots directory doesn't exit, make one.
            if not os.path.exists('./plots'): os.mkdir('./plots')
            
            #...Set up file type:
            file_type = '.pdf'
            if eps_tag == True: file_type = '.eps'
            
            #...Save pdf:
            fig.savefig('./plots/'+save_tag+'_'+self.merger_tag+file_type)
        
        elif save_frame == True:
        
            #...If saved frames directory doesn't exit, make one.
            if not os.path.exists('./frames'): os.mkdir('./frames')

            #...Set up file type:
            file_type = '.png'

            #...Save pngs:
            for ifr in np.arange(6):
                fig.savefig('./frames/frame_'+str((self.isnap)*6+ifr).zfill(4)+'.png')

        #...Report and close
        print "SAVE: ", save_tag+'_'+self.merger_tag
        plt.close()
        
    ######################################################

    def save_hdf5(self):

        #...Restrict to runs on Odyssey:
        if self.platform == 'odyssey':
            
            #...Set up hdf5 data directory path::
            hdf5_path = './hdf5_files/'+str(self.merger_tag)+'/'
        
            #...If hdf5 data directory doesn't exit, make one:
            if not os.path.exists('./hdf5_files/'): os.mkdir('./hdf5_files/')
            if not os.path.exists('./hdf5_files/'+str(self.merger_tag)+'/'): os.mkdir('./hdf5_files/'+str(self.merger_tag)+'/')

            #...Set up saving tag:
            save_tag  = self.merger_tag+'_snapshot_'+str(self.isnap).zfill(3)
            hdf5_file = hdf5_path+save_tag+'.hdf5'
            print "SAVE: ", hdf5_file

            #...Open hdf5 file:
            if self.platform == 'odyssey': data_hdf5 = h5py.File(hdf5_file,'w')
            
            #...Save gas and stellar data:
            for comp in ['gas','star']:
            
                #...Save basic data:
                data_hdf5[comp+'_id']   = np.asarray(self.all_data[comp+'_id'])
                data_hdf5[comp+'_pos']  = np.asarray(self.all_data[comp+'_pos'])
                data_hdf5[comp+'_vel']  = np.asarray(self.all_data[comp+'_vel'])
                data_hdf5[comp+'_mass'] = np.asarray(self.all_data[comp+'_mass'])
                data_hdf5[comp+'_sfr']  = np.asarray(self.all_data[comp+'_sfr'])
                #if comp == 'gas': print "saved hdf5: ", np.asarray(data_hdf5[comp+'_pos'])

                #...Save acceleration data (gas only):
                if comp == 'gas':

                    data_hdf5[comp+'_acc']       = np.asarray(self.all_data[comp+'_acc'])
                    data_hdf5[comp+'_hydro_acc'] = np.asarray(self.all_data[comp+'_hydro_acc'])
                    data_hdf5[comp+'_grav_acc']  = np.asarray(self.all_data[comp+'_grav_acc'])
                        
                    #...Save black hole data:
                    comp = 'bh'
                    data_hdf5[comp+'_id']   = np.asarray(self.all_data[comp+'_id'])
                    data_hdf5[comp+'_pos']  = np.asarray(self.all_data[comp+'_pos'])
                    data_hdf5[comp+'_vel']  = np.asarray(self.all_data[comp+'_vel'])
        
                    #...Delete data container:
                    del data_hdf5

        #...Announce error:
        elif self.platform == 'local':
            print "ERROR: You shouldn't be calling _save_hdf5() locally."
            print "ALTERNATIVE: Please sftp hdf5 files from Odyssey instead."

    ######################### PRIVATE FUNCTIONS ###########################
    
    #---Private Class Function: Sets up map.
    def _setup_fig_map(self, num=1, save_frame=False, field='mass'):
        
        #...Declare figure:
        fig = plt.figure(num, figsize=(5,5))
        
        #...Set up axes:
        ax = setup_generic_axis(fig, save_frame=save_frame)
        
        #...Set up ranges:
        ax.set_xlim([-10, 10])
        ax.set_ylim([-10, 10])
        
        #...Set up labels:
        ax.set_xlabel(r'$x$ $\rm [kpc]$', fontsize=14)
        ax.set_ylabel(r'$y$ $\rm [kpc]$', fontsize=14)
        
        return fig,ax

######################## EXTERNAL FUNCTIONS ###########################

#---Function: Generic for setting up plot axes:
def setup_generic_axis(fig, n_colors = 1, save_frame=False):

# CLIENTS: _setup_dsfr_fig()
    
    #...Declare subplot:
    ax = fig.add_subplot(1,1,1)
    
    #...Specify boundaries:
    if   save_frame == False: fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.12)
    elif save_frame ==  True: fig.subplots_adjust(left=0.0, right=1., top=1., bottom=0.)
    
    #...Specify font size of tick labes:
    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)

    #...Set up colors:
    cm        = plt.get_cmap('nipy_spectral')
    cNorm     = colors.Normalize(vmin=0, vmax=n_colors)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

    #...Set one color per curve plotted:
    ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(n_colors)])

    return ax

######################################################

#---Function: Calculates magnitude of 3xn arrays:
def calc_mag(vec):
    
    return np.sqrt(vec[:,0]*vec[:,0]+vec[:,1]*vec[:,1]+vec[:,2]*vec[:,2])

######################################################

#---Function: Calculates dot product of two 3xn arrays:
def calc_dot(vec1, vec2):

    return vec1[:,0]*vec2[:,0]+vec1[:,1]*vec2[:,1]+vec1[:,2]*vec2[:,2]

######################################################

#---Function: Calculates 2D-magnitude of 3xn arrays:
def calc_mag_xy(vec):
    
    return np.sqrt(vec[:,0]*vec[:,0]+vec[:,1]*vec[:,1])

######################################################

#---Function: Calculates 1D-magnitude of 3xn arrays:
def calc_mag_z(vec):
    
    return np.sqrt(vec[:,2]*vec[:,2])

######################################################

#---Function: Calculates cross product of two 3xn arrays:
def calc_cross(vec1, vec2):

    vec3 = np.zeros_like(vec1)
    vec3[:,0] = vec1[:,1]*vec2[:,2]-vec1[:,2]*vec2[:,1]
    vec3[:,1] = vec1[:,2]*vec2[:,0]-vec1[:,0]*vec2[:,2]
    vec3[:,2] = vec1[:,0]*vec2[:,1]-vec1[:,1]*vec2[:,0]

    return vec3

######################################################

#---Function: Calculates product of a 1xn array with a 3xn array:
def calc_prod_scal_vec(scal, vec):

    new_vec = np.zeros_like(vec)
    new_vec[:,0] = scal[:]*vec[:,0]
    new_vec[:,1] = scal[:]*vec[:,1]
    new_vec[:,2] = scal[:]*vec[:,2]
    
    return new_vec

######################################################

#---Function: Calculates unit vector from a 3xn array:
def calc_unit_vec(vec):

    mag_vec = calc_mag(vec)
    new_vec = np.zeros_like(vec)
    new_vec[:,0] = vec[:,0]/mag_vec[:]
    new_vec[:,1] = vec[:,1]/mag_vec[:]
    new_vec[:,2] = vec[:,2]/mag_vec[:]

    return new_vec

######################################################

#---Function: Rotates 3-array about z-axis:
def simple_rot_xy(vec, angle):

    new_vec = np.zeros_like(vec)
    new_vec[0] = vec[0]*np.cos(angle)-vec[1]*np.sin(angle)
    new_vec[1] = vec[0]*np.sin(angle)+vec[1]*np.cos(angle)
    new_vec[2] = vec[2]

    return new_vec

######################################################

#---Function: Rotates 3-array about y-axis:
def simple_rot_zx(vec, angle):
    
    new_vec = np.zeros_like(vec)
    new_vec[2] = vec[2]*np.cos(angle)-vec[0]*np.sin(angle)
    new_vec[0] = vec[2]*np.sin(angle)+vec[0]*np.cos(angle)
    new_vec[1] = vec[1]

    return new_vec

######################################################

#---Function: Rotates 3xn array about z-axis:
def rot_xy(vec, angle):

    new_vec = np.zeros_like(vec)
    new_vec[:,0] = vec[:,0]*np.cos(angle)-vec[:,1]*np.sin(angle)
    new_vec[:,1] = vec[:,0]*np.sin(angle)+vec[:,1]*np.cos(angle)
    new_vec[:,2] = vec[:,2]

    return new_vec

######################################################

#---Function: Rotates 3xn array about y-axis:
def rot_zx(vec, angle):
    
    new_vec = np.zeros_like(vec)
    new_vec[:,2] = vec[:,2]*np.cos(angle)-vec[:,0]*np.sin(angle)
    new_vec[:,0] = vec[:,2]*np.sin(angle)+vec[:,0]*np.cos(angle)
    new_vec[:,1] = vec[:,1]

    return new_vec

############################### END ###################################

#---Decommissioned paths (75 merger simulations):

#...Set up paths to simulations:
#            ptorrey_path = '/n/home01/ptorrey/Runs'
#            sim_paths = ['/FIRE_Mergers/Runs/Merg',
#                         '/FIRE_Mergers/Runs/Merg_lr',
#                         '/FIRE_Mergers/Runs/Merg_lr_sh03',
#                         '/FIRE_Mergers/Runs/Merg_lr_gadget',
#                         '/DavePatton/ExploratoryRuns2'
#                        ]
#
#            #...Set up simulation path:
#            self.sim_path = ptorrey_path+sim_paths[isim]
#
#            #...Set up path to orbit:
#            self.orbit_path = '/Control'
#            if self.iorbit > 0: self.orbit_path = '/Merger_'+str(self.iorbit)
#            self.orbit_path = self.sim_path+self.orbit_path
#
#            #...Set up path to snapshot
#            self.snap_path = self.orbit_path+'/output/snapshot_'+str(int(round(self.isnap/self.factor))).zfill(3)+'.hdf5'
#            print self.snap_path


