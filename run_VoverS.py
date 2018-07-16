import numpy as np
from v_over_sigma_AGS import v_over_sig
import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot as plt

irun1 = 1
irun2 = 2
irun3 = 3
irun4 = 4

snap1 = np.arange(122) #for runs with gas
index = snap1 > 29
snap2 = snap1[index] #After gas removal

tim1 = []
tim2 = []
tim3 = []
tim4 = []
Vover1 = []
Vover2 = []
Vover3 = []
Vover4 = []

tim1a = []
tim2a = []
tim3a = []
tim4a = []
Vover1a = []
Vover2a = []
Vover3a = []
Vover4a = []

##################################################  GAS RUNS
print '--------------------------------------------------'

for i in snap1:
    print 'RUN 1:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'
    
    t = i*0.05
    tim1.append(t)
    out1 = v_over_sig(irun1, i, 0)
    Vover1.append(out1[0])

for i in snap1:
    print 'RUN 2:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim2.append(t)
    out2 = v_over_sig(irun2, i, 0)
    Vover2.append(out2[0])

for i in snap1:
    print 'RUN 3:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim3.append(t)
    out3 = v_over_sig(irun3, i, 0)
    Vover3.append(out3[0])

for i in snap1:
    print 'RUN 4:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim4.append(t)
    out4 = v_over_sig(irun4, i, 0)
    Vover4.append(out4[0])

##################################################  NO GAS RUNS

for i in snap2:
    print 'RUN 1a:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim1a.append(t)
    out1a = v_over_sig(irun1, i, 1)
    Vover1a.append(out1a[0])

for i in snap2:
    print 'RUN 2a:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim2a.append(t)
    out2a = v_over_sig(irun2, i, 1)
    Vover2a.append(out2a[0])

for i in snap2:
    print 'RUN 3a:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim3a.append(t)
    out3a = v_over_sig(irun3, i, 1)
    Vover3a.append(out3a[0])

for i in snap2:
    print 'RUN 4a:'
    print ''
    print 'calculating snapshot #'+str(i)
    print '--------------------------------------------------'

    t = i*0.05
    tim4a.append(t)
    out4a = v_over_sig(irun4, i, 1)
    Vover4a.append(out4a[0])

#...For runs with gas 
ind = np.where(~np.isnan(Vover1))[0]
time1 = np.array(tim1)[ind]
VoverS1 = np.array(Vover1)[ind]

ind = np.where(~np.isnan(Vover2))[0]
time2 = np.array(tim2)[ind]
VoverS2 = np.array(Vover2)[ind]

ind = np.where(~np.isnan(Vover3))[0]
time3 = np.array(tim3)[ind]
VoverS3 = np.array(Vover3)[ind]

ind = np.where(~np.isnan(Vover4))[0]
time4 = np.array(tim4)[ind]
VoverS4 = np.array(Vover4)[ind]

#...For runs without gas
ind = np.where(~np.isnan(Vover1a))[0]
time1a = np.array(tim1a)[ind]
VoverS1a = np.array(Vover1a)[ind]

ind = np.where(~np.isnan(Vover2a))[0]
time2a = np.array(tim2a)[ind]
VoverS2a = np.array(Vover2a)[ind]

ind = np.where(~np.isnan(Vover3a))[0]
time3a = np.array(tim3a)[ind]
VoverS3a = np.array(Vover3a)[ind]

ind = np.where(~np.isnan(Vover4a))[0]
time4a = np.array(tim4a)[ind]
VoverS4a = np.array(Vover4a)[ind]



# For plotting
#plt.rcParams['lines.linewidth']=2.5#4
#plt.rcParams['axes.linewidth']=2.0#3
#plt.rcParams['lines.markersize']=5
#plt.rcParams['xtick.major.size'] = 10
#plt.rcParams['xtick.major.width'] = 3
#plt.rcParams['xtick.major.pad'] = 5
#plt.rcParams['xtick.minor.size'] = 7
#plt.rcParams['xtick.minor.width'] = 3
#plt.rcParams['xtick.minor.pad'] = 5
#plt.rcParams['ytick.major.size'] = 10
#plt.rcParams['ytick.major.width'] = 3
#plt.rcParams['ytick.major.pad'] = 5
#plt.rcParams['ytick.minor.size'] = 7
#plt.rcParams['ytick.minor.width'] = 3
#plt.rcParams['ytick.minor.pad'] = 5
#plt.rcParams['font.size']=20
#plt.rcParams['xtick.labelsize']=15
#plt.rcParams['ytick.labelsize']=15

##############################################   PLOT EVOLUTION  ########################################################

#...Set up figure box:
fig = plt.figure(figsize = (12,12))

#################################  Panel 1
ax = fig.add_subplot(111)

#...Plots V over sig vs time: WITH GAS
plt.plot(time1, VoverS1, linewidth=2.5, color = 'b')
plt.plot(time2, VoverS2, linewidth=2.5, color = 'r')
plt.plot(time3, VoverS3, linewidth=2.5, color = 'lime')
plt.plot(time4, VoverS4, linewidth=2.5, color = 'orange')

#...Plots V over sig vs time: WITHOUT GAS
plt.plot(time1a, VoverS1a, linewidth=2.5, color = 'b', linestyle = '--')
plt.plot(time2a, VoverS2a, linewidth=2.5, color = 'r', linestyle = '--')
plt.plot(time3a, VoverS3a, linewidth=2.5, color = 'lime', linestyle = '--')
plt.plot(time4a, VoverS4a, linewidth=2.5, color = 'orange', linestyle = '--')

#...Other Dwarfs
plt.axhline(y=1.7, xmin=0.0, xmax=0.3,color='lightblue') # Aquarius
plt.axhline(y=1.99, xmin=0.0, xmax=0.3,color='lightblue') # LeoA
plt.axhline(y=1.43, xmin=0.0, xmax=0.3,color='lightblue') # Pegasus
plt.axhline(y=1.01, xmin=0.0, xmax=0.3,color='lightblue') # WLM
plt.text(0.35, 1.65, r'Aquarius ', size=14)#fontdict=font)
plt.text(0.35, 1.95, r'LeoA ', size=14)
plt.text(0.35, 1.38, r'Pegaus ', size=14)
plt.text(0.35, 0.95, r'WLM ', size=14)

#...Set up axis labels:
plt.yticks(np.arange(0, 2, 0.2))
plt.xlabel(r'$time\,[Gyr]$', fontsize=20)
plt.ylabel(r'$v_{rot}/\sigma_{*}$', fontsize=20)
plt.ylim([0,2])
plt.xlim([0,6])

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2.0)

#...Legend:
#plt.legend(['medgas4', 'medgas5', 'medgas6'], loc='upper right')

#################################################  SAVING ########################################################

save_fig_file = '/data8/data/mercadf1/output/Pegasus/high_res/low_gas/run_1/Plots/lowgas_Vsig_time.png'
#save_fig_file = '/data25/rouge/gonzaa11/francisco/outputs/runmed6/Plots'
#...Report saving:
print "Saving : "
print str(save_fig_file)
    
#...Save Figure:
plt.savefig(save_fig_file)
