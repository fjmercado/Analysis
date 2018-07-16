import numpy as np
import h5py
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt
import math
from operator import truediv

#...Number of snaps
snaps = 120

#...Chose Run
irun = 4

#...Particle Type
iPartType = 4

#...Path to data
base = '/data8/data/mercadf1/output/Pegasus/high_res/run_'+ str(irun)

#...Tag for Saving
tag_list = ['Gas', 'DM', '?', '',  'Stars', '' ]
tag_type = tag_list[iPartType]

KE = np.zeros(snaps+1)
U = np.zeros(snaps+1)
L = np.zeros(snaps+1)
Lz = np.zeros(snaps+1)
time = np.zeros(snaps+1)

for i in range(snaps+1):
    #...Time
    time[i] = i*0.05
    #...Load File
    f = h5py.File(base + '/snapshot_' + str(i).zfill(3) +'.hdf5')
    
    
    #...Load Data
    xyz = np.array(f['PartType'+str(iPartType)+'/Coordinates'])

    v_xyz = np.array(f['PartType'+str(iPartType)+'/Velocities'])
    v_x = v_xyz[:,0]
    v_y = v_xyz[:,1]
    v_z = v_xyz[:,2]

    mass = np.array(f['PartType'+str(iPartType)+'/Masses'])#/(10**10)

    #...Calculate total vel
    v_tot = np.sqrt(v_x*v_x + v_y*v_y + v_z*v_z)

    #... KE
    KE[i] = 0.5*np.sum(mass*v_tot*v_tot)
    #... U                
    U[i] = -2*KE[i]
    #... L
    Ls = np.cross(xyz,v_xyz)#specific L per particle
    Ls = Ls*mass[0]
    Lst = Ls.sum(axis=0)
    Lst = Lst/sum(mass)
    L[i] = np.sqrt(Lst[0]*Lst[0] + Lst[1]*Lst[1] + Lst[2]*Lst[2])
    Lz[i] = Lst[2]

##############################################   PLOT PROFILES  ########################################################

#...Set up figure box:
fig = plt.figure(figsize = (12,12))

#################################  Panel 1

ax = fig.add_subplot(3, 1, 1)

#...Plots KE vs time:
plt.plot(time, KE, linewidth=2.5, color = 'b')#'0.75')


#...Set up axis labels:
plt.ylabel(r'$KE$', fontsize=20)
#plt.ylim([-0.01,0.01])

#...Legend:
plt.legend(['$KE$'], loc='upper right')

#################################  Panel 2

ax = fig.add_subplot(3, 1, 2)

#...Plots U vs time:
plt.plot(time, U, linewidth=2.5, color = 'k')

#...Set up axis labels:
plt.ylabel(r'$U$', fontsize=20)
#plt.ylim([-0.01,0.01])

#...Legend:
plt.legend(['$U$'], loc='upper right')

#################################  Panel 3

ax = fig.add_subplot(3, 1, 3)

#...Plots U vs time:
plt.plot(time, L, linewidth=2.5, color = 'm')
plt.plot(time, Lz, linewidth=2.5, color = 'c') 
#...Set up axis labels:
plt.xlabel(r'$time\,[Gyr]$', fontsize=20)
plt.ylabel(r'$L,\,L_z$', fontsize=20)
plt.ylim([20,30])

#...Legend:
plt.legend(['$L$', '$L_z$'], loc='upper right')

#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/high_res/run_'+ str(irun)+'/Plots/Energy_run'+str(irun)+'.png'
#save_fig_file = 'Vsig_run'+str(irun2)+'.png'
#...Report saving:
print "Saving : "
print str(save_fig_file)

#...Save Figure:
plt.savefig(save_fig_file)
