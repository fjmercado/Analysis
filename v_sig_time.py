import numpy as np
from v_over_sigma import v_sig_3D
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt

irun1 = 2
irun2 = 3

snap = np.arange(150)

tim = []
Vover1 = []
Vover2 = []

for i in snap:
    print '--------------------------------------------------'
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'
    t=(i+1.)/100.
    tim.append(t)

    out1 = v_sig_3D(irun1, i)
    print 'v_sig1: '+ str(out1)
    Vover1.append(out1)
    
    print ''

    out2 = v_sig_3D(irun2, i)
    print 'v_sig2: '+ str(out2)
    Vover2.append(out2)


ind = np.where(~np.isnan(Vover1))[0]
time1 = np.array(tim)[ind]
VoverS1 = np.array(Vover1)[ind]

ind = np.where(~np.isnan(Vover2))[0]
time2 = np.array(tim)[ind]
VoverS2 = np.array(Vover2)[ind]


##############################################   PLOT PROFILES  ########################################################

#...Set up figure box:
fig = plt.figure(figsize = (12,9))

#################################  Panel 1

ax = fig.add_subplot(1,1,1)

#...Plots V over sig vs time:
plt.plot(time1, VoverS1, linewidth=2.5, color = '0.75')
plt.plot(time2, VoverS2, linewidth=2.5, color = 'k')

#...Set up axis labels:
plt.xlabel(r'$time\,[Gyr]$', fontsize=20)
plt.ylabel(r'$v_{rot}/\sigma$', fontsize=20)
plt.ylim([0,1.0])
#...Legend:
plt.legend(['Run'+str(irun1), 'Run '+str(irun2)], loc='upper right')


#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/high_res/run_'+ str(irun2)+'/Plots/Vsig_run'+str(irun2)+'.png'
#save_fig_file = 'Vsig_run'+str(irun2)+'.png'
#...Report saving:
print "Saving : "
print str(save_fig_file)

#...Save Figure:
plt.savefig(save_fig_file)
