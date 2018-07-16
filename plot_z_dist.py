import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from pylab import *
from matplotlib import pyplot as plt

list = ['b','c','d','e','f','g','h','i','j','k','l','m']
for i in range(len(list)):
    data = np.loadtxt('10'+list[i]+'.txt')

######################  PLOT  ###################
    fig = plt.figure(figsize = (12,12))
    rc('axes',linewidth=3)
    plt.yticks(fontsize = 25)
    plt.xticks(fontsize = 25)
    plt.tick_params(which='minor',width=2,length=5)
    plt.tick_params(which='major',width=2,length=10)

    weights = np.ones_like(data)/float(len(data))
    plt.hist(data, bins = 25, weights = weights, histtype = 'step')

    plt.title('m10'+list[i], fontsize = 30)
    plt.ylabel('fraction', fontsize = 25)
    plt.xlabel('[Fe/H]', fontsize = 25)

######################  SAVE  ###################

    save_fig_file = '10'+list[i]+'.png'

    #...Report saving:
    print "Saving : "
    print str(save_fig_file)

    #...Save Figure:
    fig.savefig(save_fig_file)

