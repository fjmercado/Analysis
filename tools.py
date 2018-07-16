'''
Created:             M/07/17/17
Last Updated:        M/07/17/17
'''
#########################################################################

import numpy as np
import h5py 
import matplotlib as plt
plt.use('Agg')

#...Set up fontsize
fontsize = 14

#########################################################################

class Tools:

    #...Class initialization: Path(s) to hdf5
    def __init__(self, irun , isnap = 0, icond = 1):

        self.irun = irun

        #...path to data
        if icond == 0: #...For my ICs
            self.base = '/data8/data/mercadf1/spheric/IC_files'
            self.fname = self.base + '/New_nfwKing_lowgas_MperP400.hdf5'


        if icond == 1: #...For regular snaps
            self.base = '/data8/data/mercadf1/output/Pegasus/high_res/aruns/run_'+ str(irun) #...Path to data
            self.fname = self.base + 'a/snapshot_' + str(isnap).zfill(3) +'.hdf5' #...Generic snapshot names
    
        self.f = h5py.File(self.fname, 'r')
    
    ###########################  Initialize Private Functions  ###########################
    
    #... Private Class function: Load Data
    def load_data(self, ipart, itype):
        
        '''
        ipart = 0 (gas), 1 (DM), 2(disk), 3(bulge), 4(stars)
        '''

        '''
        itype = 0 (masses), 1 (coordinates), 2 (velocities), 3 (sfr), 4 (IDs) 
        '''

        #... Determine the particle type
        if ipart == 0:
            comp = self.f['PartType0']
        if ipart == 1:
            comp  = self.f['PartType1']
        if ipart == 2:
            comp = self.f['PartType2']
        if ipart == 3:
            comp = self.f['PartType3']
        if ipart == 4:
            comp = self.f['PartType4']

        #... Determine the data type 
        
        if itype == 0:
            dat = np.array(comp['Masses'])
        if itype == 1:
            dat = np.array(comp['Coordinates'])
        if itype == 2:
            dat = np.array(comp['Velocities'])
        if itype == 3:
            dat = np.array(comp['SFR']) #... This is probably not correct (fix it later)
        if itype == 4:
            dat = np.array(comp['IDs'])
        
        return dat

    ###########################  Initialize Public Functions  ###########################
    
    #####################################################################################
    ################################  Center of Mass  ###################################
    #####################################################################################
    
    '''
    Class Function: Determines the center of mass of distribution of particles
    '''
    def COM(self, mass, coords):
        mass_tot = sum(mass)
        mp = 0
        for j in range (0, len(mass)):
            mp = mass[j] * coords[j,:] + mp
        com = mp/mass_tot
        return com

    #####################################################################################
    #############################  Find Half-mass Radius  ###############################
    #####################################################################################

    '''
    #...Class Function: Test Private Function
    '''
    def half_mass(self, mass, coords):
        
        #...Load DM data for recentering
        DM_mass = self.load_data(1,0)
        DM_coords = self.load_data(1,1)

        #...Center of Mass Calculation (DM Center of Mass) ****MAKE A FUNCTION OF THIS???
        DM_cm = self.COM(DM_mass, DM_coords)
        
        #...Recenter:
        coords = coords - DM_cm

        #...Creates an array of particle distance from center:
        rad = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)

        #...Get rid of hubble param
        rad = rad * 0.7

        #...Calculate profile
        hist, bin_edges = np.histogram(rad, bins = 200, weights = mass, range = [0,5.0])
        hist_tot = sum(hist)*1.0
        hist = hist.cumsum()/hist_tot
        h = 0
        i = 0
        while h < 0.5:
            h = hist[i]
            half_rad = bin_edges[i]
            i = i + 1
        print 'half rad: '+str(half_rad)
        return half_rad, DM_cm

    #####################################################################################
    ##################################  Baryon Fraction  ###################################
    #####################################################################################

    '''
    Class Function: calculates the ratio of baryon mass to dynamical (total) mass within 
    the "half-light radius" (stellar half mass radius)
    '''
    def bar_frac(self):
        #...Load data
        star_coords = self.load_data(4,1)
        star_mass = self.load_data(4,0)
        gas_coords = self.load_data(0,1)
        gas_mass = self.load_data(0,0)
        DM_coords = self.load_data(1,1)
        DM_mass = self.load_data(1,0)
        
        #...Stellar half mass
        half_rad, DM_cm = self.half_mass(star_mass, star_coords)
        print 'Steller half-mass radius: '+str(half_rad)+'kpc'

        #...Resenter components
        star_coords = star_coords - DM_cm
        gas_coords = gas_coords - DM_cm

        #...Particle distances for center
        r_star = np.sqrt(star_coords[:,0]**2 + star_coords[:,1]**2 + star_coords[:,2]**2)*0.7
        r_gas = np.sqrt(gas_coords[:,0]**2 + gas_coords[:,1]**2 + gas_coords[:,2]**2)*0.7
        r_DM = np.sqrt(DM_coords[:,0]**2 + DM_coords[:,1]**2 + DM_coords[:,2]**2)*0.7

        #...Create masks
        mask_s = np.where(r_star <= half_rad)[0]
        mask_g = np.where(r_gas <= half_rad)[0]
        mask_d = np.where(r_DM <= half_rad)[0]

        #...New masses
        star_mass = star_mass[mask_s]
        gas_mass = gas_mass[mask_g]
        DM_mass = DM_mass[mask_d]

        #...Sum of component masses
        print 'Within the stellar half mass radius:'
        smass = sum(star_mass)*10**10
        print '    Stellar mass = '+str(smass)+' Msun'
        gmass = sum(gas_mass)*10**10
        print '    Gas mass = '+str(gmass)+' Msun'
        dmass = sum(DM_mass)*10**10
        print '    DM mass = '+str(dmass)+' Msun'

        #...Gas Fracstion
        frac = (gmass+smass)/(smass+gmass+dmass)
        print ''
        print 'Baryon fraction: '+str(frac)
        return frac
