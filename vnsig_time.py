import numpy as np
from v_over_sigma_AGS import v_over_sig
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt

irun = 4

snap = np.arange(120)

tim = []
v_rot = []
sig = []

print '--------------------------------------------------'

for i in snap:
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim.append(t)
    out1 = v_over_sig(irun, i)
    v_rot.append(out1[2])
    sig.append(out1[3])

ind = np.where(~np.isnan(v_rot))[0]
time1 = np.array(tim)[ind]
vrot = np.array(v_rot)[ind]

ind = np.where(~np.isnan(sig))[0]
time2 = np.array(tim)[ind]
sigma = np.array(sig)[ind]

##############################################   PLOT PROFILES  ########################################################

#...Set up figure box:
fig = plt.figure(figsize = (12,12))

#################################  Panel 1

#ax = fig.add_subplot(2, 1, 1)

#...Plots V over sig vs time:
plt.plot(time1, vrot, linewidth=2.5, color = 'r')#'0.75')
plt.plot(time2, sigma, linewidth=2.5, color = 'k')

#...Set up axis labels:
plt.xlabel(r'$time\,[Gyr]$', fontsize=20)
plt.ylabel(r'$v_{rot},\,sigma\,[km\,s^{-1}]$', fontsize=20)
#plt.ylim([0,50.0])
#...Legend:
plt.legend(['$V_{rot}$', '$Sigma$'], loc='upper right')


#################################  Panel 2
#ax = fig.add_subplot(2, 1, 2)
#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/high_res/run_'+ str(irun)+'/Plots/VnS_run'+str(irun)+'.png'
#save_fig_file = 'Vsig_run'+str(irun2)+'.png'
#...Report saving:
print "Saving : "
print str(save_fig_file)

#...Save Figure:
plt.savefig(save_fig_file)
