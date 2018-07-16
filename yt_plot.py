import yt
import yt.units as units
import pylab
from yt.analysis_modules.level_sets.api import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

"*************** LOAD DATA ********************"
#...Choose snapshot:
isnap = 100
 
#...Chose run:
run = 'run_3' 
 
# Name and location of Snapshot (data)
fname = '/data8/data/mercadf1/output/Pegasus/'+run+'/snapshot_' + str(isnap).zfill(3) +'.hdf5'
 
# Width of whatever plot
my_width = 15.0

# Color map you'd like to use
my_cmap = 'GRN-RED-BLU-WHT'
 
# The limit of the bounding box being created
bbox_lim = 1e5 #kpc


# Define center 
halocenter = (0,0,0)
 
# Create the bounding box that will "house" the data
bbox = [[-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim],
        [-bbox_lim,bbox_lim]]
 
# Loading the data within the bounding box
ds = yt.load(fname, bounding_box=bbox)
 
# Load the data from the Yt Sample Data
#ds = yt.load("/Users/student/Desktop/IsolatedGalaxy/galaxy0030/galaxy0030")
ds.index
ad= ds.all_data()


#define geometry
sp = ds.sphere(halocenter, (bbox_lim, 'kpc'))
"*********** CREATE AND SAVE SLICE PLOTS ******************"

## Plot of an X-axis Slice
#slc = yt.SlicePlot(ds, 'x', [('gas', 'density')], width = (my_width, 'kpc'))
#slc.set_cmap(field="density", cmap= my_cmap)
#slc.set_zlim('density', 10**-28, 10**-23)
#slc.set_ylim(-100,100)
#slc.set_xlim(-100,100)
#slc.save('/data8/data/mercadf1/output/Pegasus/'+run+'/Plots/GasMassSlice_x_snapshot_' + str(isnap).zfill(3) +'.png')

#
## Plot of an Y-axis Slice
#slc = yt.SlicePlot(ds, 'y', [('gas', 'density')], width = (my_width, 'kpc'))
#slc.set_cmap(field="density", cmap= my_cmap)
#slc.save('/Users/francisco/DATA/Pegasus/'+run+'/plots/GasDensitySlice_y_snapshot_' + str(isnap).zfill(3) +'.png')
#
## Plot of a Z-axis Slice
#slc = yt.SlicePlot(ds, 'z', [('gas', 'density')], width = (my_width, 'kpc'))
#slc.set_cmap(field="density", cmap= my_cmap)
#slc.save('/Users/francisco/DATA/Pegasus/'+run+'/plots/GasDensitySlice_z_snapshot_' + str(isnap).zfill(3) +'.png')
 


"************ CREATE AND SAVE A PARTICLE PLOT ***************"

# create our plot
#p = yt.ParticleProjectionPlot(ds,2,[('PartType1','Masses')], width = (5,5), center = halocenter, data_source = sp)
#p.set_cmap(field="particle_mass", cmap= my_cmap)
#p.set_unit('particle_mass', 'Msun')
#save result
#p.save('/Users/francisco/DATA/Pegasus/run_2/plots/particle_plot_snapshot_' + str(isnap).zfill(3) +'.png')

#####...Test Particle Plot
p0 = yt.ParticlePlot(ds, ('PartType0','particle_position_x'), ('PartType0','particle_position_y'),('PartType0', 'density'))
#p0.set_zlim('density', 5*10**-6, 10**-3 )
p0.save('/data8/data/mercadf1/output/Pegasus/'+run+'/Plots/Gas_ParticlePlot_x_snapshot_' + str(isnap).zfill(3) +'.png')

#p1 = yt.ParticlePlot(ds, ('PartType1','particle_position_x'), ('PartType1','particle_position_y'),('PartType1', 'particle_mass'))
#p1.set_zlim('particle_mass', 5*10**-6, 10**-3 )
#p1.save('/data8/data/mercadf1/output/Pegasus/'+run+'/Plots/DM_ParticlePlot_x_snapshot_' + str(isnap).zfill(3) +'.png')

#p2 = yt.ParticlePlot(ds, ('PartType4','particle_position_x'), ('PartType4','particle_position_y'),('PartType4', 'particle_mass'))
#p2.set_zlim('particle_mass', 5*10**-6, 10**-3 )
#p2.save('/data8/data/mercadf1/output/Pegasus/'+run+'/Plots/Stars_ParticlePlot_x_snapshot_' + str(isnap).zfill(3) +'.png')

#p1 = yt.ParticlePlot(ds,[('ParticleType4','particle_position_x'),('ParticleType4','particle_position_y')],('ParticleType4','particle_mass'), width = (my_width, 'kpc'))
#p1.set_cmap(field="particle_mass", cmap= my_cmap)
#p1.set_unit('particle_mass', 'Msun')
##save result
#p1.save('/Users/francisco/DATA/Pegasus/run_2/plots/particle_plot_snapshot_' + str(isnap).zfill(3) +'.png')

"************ CREATE OFF-AXIS PROJECTION PLOT *****************"
 
## Create a 15 kpc radius sphere, centered on the center of the sim volume
#sp = ds.sphere("center", (my_width, "kpc"))
#
## Get the angular momentum vector for the sphere.
#L = sp.quantities.angular_momentum_vector()
#
#print("Angular momentum vector: {0}".format(L))
#
## Create an OffAxisProjectionPlot of density centered on the object with the L
## vector as its normal and a width of 25 kpc on a side
#p = yt.OffAxisProjectionPlot(ds, L, "density", sp.center, (25, "kpc"))
#p.set_cmap(field="density", cmap= my_cmap)
#p.save('/Users/francisco/DATA/Pegasus/run_1/plots/OffAxizProjPlot_snapshot_' + str(isnap).zfill(3) +'.png')
