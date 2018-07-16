import h5py
import numpy as np
import math as m
###############################################  BOOKKEEPING  #####################################################

#...Choose snapshot
#isnap = 9

#...Choose run:
#irun = 6

#############################################  LOAD DATA  #########################################################


#def v_over_sig(dat,isnap, plane='xy', nbins = 30, pt = 'Star', X = 3.0):
def v_over_sig(irun, isnap, run_Type, plane='xy', nbins = 30, pt = 'Star', ex = 3.0):
    
    #irun = 6

    #...Path to data
    if run_Type == 0:
        base = '/data8/data/mercadf1/output/Pegasus/high_res/low_gas/run_'+ str(irun)
    #base = '/data25/rouge/gonzaa11/francisco/outputs/runmed'+ str(irun)
    if run_Type == 1:
        base = '/data8/data/mercadf1/output/Pegasus/high_res/low_gas/run_'+ str(irun)+'a'
#...Load File
    f = h5py.File(base + '/snapshot_' + str(isnap).zfill(3) +'.hdf5')

    p = np.array(f['PartType4/Coordinates'])
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]

    v = np.array(f['PartType4/Velocities'])
    v_x = v[:,0]
    v_y = v[:,1]
    v_z = v[:,2]

    ################### Cylindrical Coordinate Stuff

    r = (x*x + y*y)**0.5
    yoverx = y/x
    phi = np.zeros(len(yoverx))
    v_r = np.zeros(len(yoverx))
    v_phi = np.zeros(len(yoverx))
    for i in range(len(phi)):
        phi[i] = np.arctan2(y[i],x[i])
        v_r[i] = v_x[i]*m.cos(phi[i]) + v_y[i]*m.sin(phi[i])
        v_phi[i] = -v_x[i]*m.sin(phi[i]) + v_y[i]*m.cos(phi[i])

    ##############
    N = len(x)
    Tratio = np.ndarray(nbins)
    q = np.ndarray(nbins)
           
   # alltheta = np.linspace(0,np.pi,nbins)
   # for i in range(nbins):
   #     rot_theta = alltheta[i]
   #     # Now get los velocity for given rotation angle and plane
   #     if plane == 'xy':
   #         vlos = v[:,2]
   #         wd = p[:,0] * np.cos(rot_theta) - p[:,1] * np.sin(rot_theta)
   #         ht = p[:,1] * np.cos(rot_theta) + p[:,0] * np.sin(rot_theta)
   #     elif plane == 'xz':
   #         vlos = -v[:,1]
   #         wd = p[:,0] * np.cos(rot_theta) - p[:,2] * np.sin(rot_theta)
   #         ht = p[:,2] *  np.cos(rot_theta) + p[:,0] * np.sin(rot_theta)
   #     elif plane == 'yz':
   #         vlos = v[:,0]
   #         wd = p[:,1] * np.cos(rot_theta) - p[:,2] * np.sin(rot_theta)
   #         ht = p[:,2] * np.cos(rot_theta) + p[:,1] * np.sin(rot_theta)
    
    
    wbins = np.linspace(0,ex,nbins+1)
    hbins = np.linspace(0,ex,nbins+1)
    rbins = np.linspace(0,ex,nbins+1)
   
    vrotb = np.ndarray(nbins)
    sig3Db = np.ndarray(nbins)
    sig1D = np.ndarray(nbins)
    ninbin = np.ndarray(nbins)
    ratiob = np.ndarray(nbins)
    rmidbin = np.ndarray(nbins)
    xmidbin = np.ndarray(nbins)
    ymidbin = np.ndarray(nbins)
    index = []
    
    for i in range(nbins):
        indbin = np.where((r >= rbins[i])&(r < rbins[i+1]))[0]
        index.append(indbin)
        ninbin[i] = len(indbin)
        rmidbin[i] = (rbins[i] + rbins[i+1]) / 2.0
        xmidbin[i] = (x[i] + x[i+1]) / 2.0
        ymidbin[i] = (y[i] + y[i+1]) / 2.0
        if ninbin[i] > 5:
            vrotb[i] = np.mean(v_phi[indbin])
            sig3Db[i] = np.sqrt((v_r[indbin].std())**2. + (v_phi[indbin].std())**2. + (v_z[indbin].std())**2.)
            sig1D[i] = sig3Db[i]/np.sqrt(3.)
            if sig3Db[i] > 0:
                ratiob[i] = vrotb[i] / sig1D[i]
            else:
                ratiob[i] = vrotb[i]
        else:
            vrotb[i] = sig1D[i] = ratiob[i] = 0
    
    Tratio, q, v_rot, sigma  = get_sauron_ratio_q(ninbin, vrotb, sig1D, ratiob, xmidbin, ymidbin)
    

    #return Tratio[np.argmin(q)], 1 - np.amin(q), alltheta[np.argmin(q)]#, N
    return Tratio, 1 - q, v_rot, sigma #,Ratio3D, N
    print index





def get_sauron_ratio_q(ninbin, vrot, sig, ratio, xmidbin, ymidbin):
    """
    Calculatea vrot / sigma by summing bin quantities
    See Cappellari et al. 2007 (SAURON)    
    
    
    """
    
    v_rot = np.mean(vrot)
    sigma = np.mean(sig)
    Tratio = np.sqrt(np.sum(ninbin * vrot**2)/np.sum(ninbin * sig**2))
    q = np.sqrt(np.sum(ninbin * ymidbin**2)/np.sum(ninbin * xmidbin**2))
    print 'v_rot: '+str(v_rot*.001)+'km/sec'
    print 'sigma: '+str(sigma)+'km/sec'
    print 'v/sigma: '+str(Tratio)
    return Tratio, q, v_rot, sigma


##### hasta aqui he revisado! Luego le sigo-----------
### creo que lo de abajo sirve para solo incluir la tansformacion a cilindricas en la parte superior


def v_sig_3D(irun, isnap):
    base = '/data8/data/mercadf1/output/Pegasus/high_res/run_'+ str(irun)

    #...Load File
    f = h5py.File(base + '/snapshot_' + str(isnap).zfill(3) +'.hdf5')

    p = np.array(f['PartType4/Coordinates'])
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]
    
    v = np.array(f['PartType4/Velocities'])
    v_x = v[:,0]
    v_y = v[:,1]
    v_z = v[:,2]
    v_tot = (v_x*v_x + v_y*v_y + v_z*v_z)**0.5
    
    #...V_phi average (cylindrical Coords)
    #...Note: For this to work, the stelar component's total angular momentum vector should already by oriented in the z direction

    ################### Cylindrical Coordinate Stuff

    r = (x*x + y*y)**0.5
    yoverx = y/x
    phi = np.zeros(len(yoverx))
    for i in range(len(phi)):
        phi[i] = m.atan(yoverx[i])

    #v_phi = -1.0 * (r*(v_x/x - v_y/y))/ (x/y + y/x)
    v_phi = -1.0*v_x*m.sin(phi)+ v_y*m.cos(phi)
    
    v_phi_avg = np.mean(v_phi) #....We use this for v/sig

    #...Sigma (std dev of velocities)
    sigma_3D = np.std(v_tot)
    
    sigma = sigma_3D/(3 ** .5)
    
    #...Calculate v/sig

    v_sig = v_phi_avg/sigma

    
    return v_sig










