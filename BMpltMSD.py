import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import time


def msdcalculator(path='/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/ecc0'):
    #return msd, msdx, msdy. i.e. msd[0]: Mean of the msd. msd[1]: SEM of the msd
    # n = 100000
    n = 3000 # number of frame for each clip
    # n_test = 5
    n_test = 50 # number of clips
    pos = np.loadtxt(path + '/pos.txt')
    # pos = np.loadtxt(path + '/pos.txt')

    msd = []
    msd_x = []
    msd_y = []
    for i in range(n_test):
        t0 = time.time()
        x = pos[i*n:(i+1)*n, 0]
        y = pos[i*n:(i+1)*n, 1]
        deno = np.arange(1, n + 1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        msd.append((2 * (np.mean(x ** 2 + y ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                                                  deno) - 2 * np.divide(
            np.correlate(y, y, mode='full'), deno)))
        msd_x.append((2 * (np.mean(x ** 2 )) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                                                  deno)))
        msd_y.append((2 * (np.mean(y ** 2)) - 2 * np.divide(np.correlate(y, y, mode='full'),
                                                            deno)))
        # print('One s')
    msd_mean = np.mean(msd, axis=0)
    msd_err = np.std(msd,axis=0)/np.sqrt(n_test)
    msd= []
    msd.append(msd_mean)
    msd.append(msd_err)
    msdx_mean = np.mean(msd_x, axis=0)
    msdx_err = np.std(msd_x,axis=0)/np.sqrt(n_test)
    msdx = []
    msdx.append(msdx_mean)
    msdx.append(msdx_err)
    msdy_mean = np.mean(msd_y, axis=0)
    msdy_err = np.std(msd_y, axis=0)/np.sqrt(n_test)
    msdy = []
    msdy.append(msdy_mean)
    msdy.append(msdy_err)
    return msd, msdx, msdy
def msd_individualclips(path):
    # return msd, msdx, msdy. i.e. msd[0]: Mean of the msd. msd[1]: SEM of the msd
    # n = 100000
    n = 3000  # number of frame for each clip
    # n_test = 5
    n_test = 50  # number of clips for clips with T4
    # n_test = 10000 # number of clips for clips without T4
    pos = np.loadtxt(path + '/pos.txt')
    # pos = np.loadtxt(path + '/pos.txt')

    msd = []
    msd_x = []
    msd_y = []
    for i in range(n_test):
        x = pos[i * n:(i + 1) * n, 0]
        y = pos[i * n:(i + 1) * n, 1]
        deno = np.arange(1, n + 1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        msd.append((2 * (np.mean(x ** 2 + y ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                                   deno) - 2 * np.divide(
            np.correlate(y, y, mode='full'), deno)))
        msd_x.append((2 * (np.mean(x ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                            deno)))
        msd_y.append((2 * (np.mean(y ** 2)) - 2 * np.divide(np.correlate(y, y, mode='full'),
                                                            deno)))

    # read the 2nd batch msd
    n = 3000  # number of frame for each clip
    # n_test = 5
    n_test = 200  # number of clips
    pos = np.loadtxt(path + '/pos2.txt')
    # pos = np.loadtxt(path + '/pos.txt')
    for i in range(n_test):
        x = pos[i * n:(i + 1) * n, 0]
        y = pos[i * n:(i + 1) * n, 1]
        deno = np.arange(1, n + 1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        msd.append((2 * (np.mean(x ** 2 + y ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                                   deno) - 2 * np.divide(
            np.correlate(y, y, mode='full'), deno)))
        msd_x.append((2 * (np.mean(x ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                            deno)))
        msd_y.append((2 * (np.mean(y ** 2)) - 2 * np.divide(np.correlate(y, y, mode='full'),
                                                            deno)))

    return msd, msd_x, msd_y
def msd_avg(item):
    ##################### Take average of the clips #############
    # This is for the calibration of simulation time scale
    # item: 'ecc0' for example
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/'+item)
    msd = np.loadtxt('msd_sim')
    msdx = np.loadtxt('msdx_sim')
    msdy = np.loadtxt('msdy_sim')
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/msdsim_calibration/Simulation/'+item)
    msd_t = np.mean(msd, axis=0)
    std_t = np.std(msd, axis=0)
    np.savetxt('msdsim_mean.txt', msd_t)
    np.savetxt('stdsim_mean.txt', std_t)

    msd_t = np.mean(msdx, axis=0)
    std_t = np.std(msdx, axis=0)
    np.savetxt('msdxsim_mean.txt', msd_t)
    np.savetxt('stdxsim_mean.txt', std_t)

    msd_t = np.mean(msdy, axis=0)
    std_t = np.std(msdy, axis=0)
    np.savetxt('msdysim_mean.txt', msd_t)
    np.savetxt('stdysim_mean.txt', std_t)
    return
def msd_individualclips_noT(path):
    # return msd, msdx, msdy. i.e. msd[0]: Mean of the msd. msd[1]: SEM of the msd
    # n = 100000
    n = 3000  # number of frame for each clip
    # n_test = 5
    # n_test = 50  # number of clips for clips with T4
    n_test = 3000 # number of clips for clips without T4
    print('Loading files...')
    pos = np.loadtxt(path + '/pos.txt', max_rows=int(n*n_test))
    print('Calculating msd...')
    # pos = np.loadtxt(path + '/pos.txt')

    msd = []
    msd_x = []
    msd_y = []
    for i in range(n_test):
        x = pos[i * n:(i + 1) * n, 0]
        y = pos[i * n:(i + 1) * n, 1]
        deno = np.arange(1, n + 1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        msd.append((2 * (np.mean(x ** 2 + y ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                                   deno) - 2 * np.divide(
            np.correlate(y, y, mode='full'), deno)))
        msd_x.append((2 * (np.mean(x ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                            deno)))
        msd_y.append((2 * (np.mean(y ** 2)) - 2 * np.divide(np.correlate(y, y, mode='full'),
                                                            deno)))
    os.chdir(path)
    np.savetxt('msd', msd)
    np.savetxt('msd_x', msd_x)
    np.savetxt('msd_y', msd_y)
    return msd, msd_x, msd_y
def msd_individualclips_ecc0newsim(path):
    # return msd, msdx, msdy. i.e. msd[0]: Mean of the msd. msd[1]: SEM of the msd
    # n = 100000
    n = 5000  # number of frame for each clip
    # n_test = 5
    n_test = 250  # number of clips for clips with T4
    # n_test = 10000 # number of clips for clips without T4
    pos = np.loadtxt(path + '/pos3.txt')
    # pos = np.loadtxt(path + '/pos.txt')

    msd = []
    msd_x = []
    msd_y = []
    print('Calculating MSD...')
    for i in range(n_test):
        x = pos[i * n:(i + 1) * n, 0]
        y = pos[i * n:(i + 1) * n, 1]
        deno = np.arange(1, n + 1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        msd.append((2 * (np.mean(x ** 2 + y ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                                   deno) - 2 * np.divide(
            np.correlate(y, y, mode='full'), deno)))
        msd_x.append((2 * (np.mean(x ** 2)) - 2 * np.divide(np.correlate(x, x, mode='full'),
                                                            deno)))
        msd_y.append((2 * (np.mean(y ** 2)) - 2 * np.divide(np.correlate(y, y, mode='full'),
                                                            deno)))

    return msd, msd_x, msd_y
# path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Calibrate/'
# path2 = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/calibration/test/'
# # msd_individualclips_noT(path)
#
# path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/' # position file directory
# os.chdir(path)
# dirlist = ['ecc0', 'ecc06','ecc08', 'ecc09', 'ecc095'] # list directory
# # dirlist = ['ecc0995']
# n = 3000 # step number
# n_test = 50 # number of simulated videos

#################### exponent fitting ##################
def exponentfit():
    ftime = 20 # fitting frames
    def fitfunc(x,a,b):# fitting function for MSD. Exponent
        return a*x**b
    sol = []
    cov = []
    msd_avg = []
    for folder in dirlist:
        msd, msdx, msdy = msdcalculator(path + folder)
        msd_avg.append(msdy[0]) # adding msd/msdx or msdy to proccessing data.
    for item in msd_avg:
        msd_tmp = item[n-1:(n-1+ftime)]
        t_tmp = np.arange(0,ftime)
        popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp, p0=[1,0.5])
        sol.append(popt)
        cov.append(pcov)
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/')
    print(sol)
    print(cov)
    f = open('msdfit.txt', 'a')
    np.savetxt(f,sol)
    f.close()
    f = open('msdfit_cov.txt', 'a')
    for item in cov:
        np.savetxt(f,item)
    f.close()
def exponentfit_individualclips():
    def fitfunc(x, a, b):  # fitting function for MSD. Exponent
        return a * x ** b
    sol = []
    cov = []
    solx = []
    covx = []
    soly = []
    covy = []
    msd_avg = []
    path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/'
    ind = 0
    for folder in datalist:
        if folder == 'ecc0':
            msd, msdx, msdy = msd_individualclips_ecc0newsim(path + folder)
        else:
            msd, msdx, msdy = msd_individualclips(path + folder)

        ftime = int(1 / ratio[ind]) # fitting frames
        ff = 0  # first frame. Corresponds to the first datapoint in experiment
        # fig = plt.figure()
        for itemidx in range(len(msd)):
            # msd
            msd_tmp = msd[itemidx]
            start = int((len(msd_tmp)-1)/2)+ff # start of the msd. Only right hand side is taken because of symmetry.
            end = int(start + ftime)
            t_tmp = np.arange(ff, ff+ftime)
            msd_tmp[start] = 1e-15
            popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp[start:end], p0=[1, 0.5])
            # if itemidx<=15:
            #     ax = fig.add_subplot(4,4,itemidx+1)
            #     ax.plot(t_tmp, msd_tmp[start:end], '+')
            #     ax.plot(t_tmp, fitfunc(t_tmp, popt[0], popt[1]), '+')
            sol.append(popt)
            cov.append(pcov)
            # msdx
            msd_tmp = msdx[itemidx]
            start = int((len(msd_tmp) - 1) / 2)  # start of the msd. Only right hand side is taken because of symmetry.
            end = int(start + ftime)
            msd_tmp[start] = 1e-15
            t_tmp = np.arange(ff, ff+ftime)
            popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp[start:end], p0=[1, 0.5])
            solx.append(popt)
            covx.append(pcov)
            # msdy
            msd_tmp = msdy[itemidx]
            msd_tmp[start] = 1e-15
            start = int((len(msd_tmp) - 1) / 2)  # start of the msd. Only right hand side is taken because of symmetry.
            end = int(start + ftime)
            t_tmp = np.arange(ff, ff+ftime)
            popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp[start:end], p0=[1, 0.5])
            soly.append(popt)
            covy.append(pcov)
        ind+=1
        plt.show()

        os.chdir(
            '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/MSDfit/CalibratedFitting/')
        # save msd
        f = open('msdfit_clip_'+folder+'.txt', 'w')
        np.savetxt(f, sol)
        f.close()
        f = open('msdfit_clip_'+folder+'_cov.txt', 'w')
        for item in cov:
            np.savetxt(f, item)
        f.close()
        # save msdx
        f = open('msdfitx_clip_'+folder+'.txt', 'w')
        np.savetxt(f, solx)
        f.close()
        f = open('msdfitx_clip_'+folder+'_cov.txt', 'w')
        for item in covx:
            np.savetxt(f, item)
        f.close()
        # save msdy
        f = open('msdfity_clip_'+folder+'.txt', 'w')
        np.savetxt(f, soly)
        f.close()
        f = open('msdfity_clip_'+folder+'_cov.txt', 'w')
        for item in covy:
            np.savetxt(f, item)
        f.close()
    return
############ Plot #########################
def plotter():
    ########### fitted lines ###########
    t_t = np.linspace(1e-2,100,100)
    param = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msd/msdfit.txt')
    xparam = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdx/msdfit.txt')
    yparam = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdy/msdfit.txt')
    # param = np.loadtxt(
    #     '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/test/msdfit_clip_ecc0.txt')
    # xparam = np.loadtxt(
    #     '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/test/msdfitx_clip_ecc0.txt')
    # yparam = np.loadtxt(
    #     '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/test/msdfity_clip_ecc0.txt')
    print('Mean exponent:' + str(np.mean(param[:, 1])) + '+-' + str(np.std(param[:, 1])))
    print('X Mean exponent:' + str(np.mean(xparam[:, 1])) + '+-' + str(np.std(xparam[:, 1])))
    print('Y Mean exponent:' + str(np.mean(yparam[:, 1])) + '+-' + str(np.std(yparam[:, 1])))
    y1 = 2e-2*(t_t/t_t[1])**0.85
    y2 = 2e-2*(t_t/t_t[1])**0.92
    y31 = 1e-2*(t_t/t_t[1])**0.95
    y32 = 4e-3*(t_t)**0.57


    fig =plt.figure(figsize=[12,10],tight_layout=True)
    ax = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    t = np.arange(1,n)
    ind = 1
    cmap = cm.get_cmap("inferno").colors
    colorint = 30
    colorshift = 0.8
    marker = ['o','v','^','<','>','s','+','x']
    for folder in dirlist:
        msd, msdx, msdy = msdcalculator(path+folder)
        ax.plot(t, msd[0][n:], color = cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=dirlist[ind-1])
        ax2.plot(t, msdx[0][n:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=dirlist[ind-1])
        ax3.plot(t, msdy[0][n:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=dirlist[ind-1])
        ind+=1
    ax.plot(t_t,y1,'--',color='r', label=r'$\alpha = 0.85$')
    ax2.plot(t_t,y2,'--',color='r', label=r'$\alpha = 0.92$')
    ax3.plot(t_t,y31,'--',color='r', label=r'$\alpha = 0.95$')
    # ax3.plot(t_t,y32,ls='dashdot',color='b', label=r'$\alpha = 0.57$')

    # t = np.arange(0,500)
    # msd_scale = 2e-2*(t)**(0.8)
    #
    # ax.errorbar(x=np.arange(1,n), y=msd_mean[int(n):], yerr=msd_err[int(n):])
    # ax.plot(t,msd_scale,'--')
    ax.set_xlim([8e-1, n+1000])
    ax.set_ylim([9e-3, 3])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'MSD ($\mu$m)',fontsize=15)
    ax.set_xlabel('Lag time',fontsize=15)
    ax.legend(fontsize=10)
    #
    #
    # ax2.errorbar(x=np.arange(1,n), y=msdx_mean[int(n):], yerr=msdx_err[int(n):])
    # ax2.plot(t,msd_scale,'--')
    ax2.set_xlim([8e-1, n+1000])
    ax2.set_ylim([9e-3, 3])
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'MSD$_{||}$ ($\mu$m)',fontsize=15)
    ax2.set_xlabel('Lag time',fontsize=15)
    ax2.legend(fontsize=10)

    #
    # ax3.errorbar(x=np.arange(1,n), y=msdy_mean[int(n):], yerr=msdy_err[int(n):])
    # ax3.plot(t,msd_scale,'--')
    ax3.set_xlim([8e-1, n+1000])
    ax3.set_ylim([4e-3, 3])
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'MSD$_{\perp}$ ($\mu$m)',fontsize=15)
    ax3.set_xlabel('Lag time',fontsize=15)
    ax3.legend(fontsize=10)

    x = np.arange(0,len(param[:,1]))
    ax4.plot(x[:-1], param[:-1,1], 's', label='MSD', markersize=10)
    ax4.plot(x[:-1], xparam[:-1,1],'o', label='MSD major', markersize=10)
    ax4.plot(x[:-1], yparam[:-1,1],'^', label='MSD minor', markersize=10)
    ax4.set_xticklabels(['', '0', '0.6', '0.8', '0.9', '0.95', '0.98', '0.995'])
    ax4.legend(loc = 'lower left', fontsize=10)
    ax4.set_ylabel(r'$\alpha$', fontsize=15)
    ax4.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
    return
def plotter_2batchdata():
    ########### fitted lines ###########
    t_t = np.linspace(1e-2,100,100)
    param = [0.69, 0.73, 0.73, 0.72, 0.73]  # Time scale calibrated fitting
    xparam = [0.68, 0.75, 0.78, 0.78, 0.81]
    yparam = [0.69, 0.71, 0.67, 0.61, 0.60]
    alpha_sem_sim = [6e-3, 3e-3,3e-3,3e-3,2e-3]
    alphax_sem_sim = [6e-3, 5e-3, 4e-3, 4e-3,
                      3e-3]  # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
    alphay_sem_sim = [6e-3, 4e-3, 4e-3, 4e-3, 4e-3]
    # print('Mean exponent:' + str(np.mean(param[:, 1])) + '+-' + str(np.std(param[:, 1])))
    # print('X Mean exponent:' + str(np.mean(xparam[:, 1])) + '+-' + str(np.std(xparam[:, 1])))
    # print('Y Mean exponent:' + str(np.mean(yparam[:, 1])) + '+-' + str(np.std(yparam[:, 1])))
    # y1 = 2e-2*(t_t/t_t[1])**0.85
    # y2 = 2e-2*(t_t/t_t[1])**0.92
    # y31 = 1e-2*(t_t/t_t[1])**0.95
    # y32 = 4e-3*(t_t)**0.57


    fig =plt.figure(figsize=[12,10],tight_layout=True)
    ax = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    ind = 1
    cmap = cm.get_cmap("inferno").colors
    colorint = 30
    colorshift = 0.8
    xlimr = 2000
    marker = ['o','v','^','<','>','s','+','x']


    for folder in datalist:
        msd, msdx, msdy = msd_individualclips(path+folder)
        n = int((len(np.mean(msd, axis=0))-1)/2)
        a = np.mean(msd,axis=0)[n:]
        t = np.arange(0, len(a))
        ax.plot(t, np.mean(msd,axis=0)[n:], color = cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=datalist[ind-1])
        ax2.plot(t, np.mean(msdx,axis=0)[n:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=datalist[ind-1])
        ax3.plot(t, np.mean(msdy,axis=0)[n:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=datalist[ind-1])

        # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t', t)
        # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd', np.mean(msd,axis=0)[n:])
        # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msdx', np.mean(msdx,axis=0)[n:])
        # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msdy',
        #            np.mean(msdy, axis=0)[n:])
        # print('a')
        ind+=1
    # ax.plot(t_t,y1,'--',color='r', label=r'$\alpha = 0.85$')
    # ax2.plot(t_t,y2,'--',color='r', label=r'$\alpha = 0.92$')
    # ax3.plot(t_t,y31,'--',color='r', label=r'$\alpha = 0.95$')
    # ax3.plot(t_t,y32,ls='dashdot',color='b', label=r'$\alpha = 0.57$')

    # t = np.arange(0,500)
    # msd_scale = 2e-2*(t)**(0.8)
    #
    # ax.errorbar(x=np.arange(1,n), y=msd_mean[int(n):], yerr=msd_err[int(n):])
    # ax.plot(t,msd_scale,'--')
    ax.set_xlim([8e-1, xlimr])
    ax.set_ylim([9e-3, 3])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$',fontsize=15)
    ax.set_xlabel('Lag time',fontsize=15)
    ax.legend(fontsize=10)
    #
    #
    # ax2.errorbar(x=np.arange(1,n), y=msdx_mean[int(n):], yerr=msdx_err[int(n):])
    # ax2.plot(t,msd_scale,'--')
    ax2.set_xlim([8e-1, xlimr])
    ax2.set_ylim([9e-3, 3])
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$',fontsize=15)
    ax2.set_xlabel('Lag time',fontsize=15)
    ax2.legend(fontsize=10)

    #
    # ax3.errorbar(x=np.arange(1,n), y=msdy_mean[int(n):], yerr=msdy_err[int(n):])
    # ax3.plot(t,msd_scale,'--')
    ax3.set_xlim([8e-1, xlimr])
    ax3.set_ylim([4e-3, 3])
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$',fontsize=15)
    ax3.set_xlabel('Lag time',fontsize=15)
    ax3.legend(fontsize=10)

    x = np.arange(0,len(param))
    ax4.errorbar(x, param, yerr = alpha_sem_sim, marker='s', label=r'$\langle(\Delta r)^2\rangle$', markersize=10,capsize=4, ls='None')
    ax4.errorbar(x, xparam,yerr = alphax_sem_sim, marker='o', label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=10, capsize=4, ls='None')
    ax4.errorbar(x, yparam,yerr = alphay_sem_sim, marker='^', label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=10, capsize=4, ls='None')

    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/param', param)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/paramsem', alpha_sem_sim)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xparam', xparam)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xparamsem', alphax_sem_sim)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/yparam', yparam)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/yparamsem', alphay_sem_sim)
    ax4.set_xticks(x)
    ax4.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax4.legend(loc = 'lower left', fontsize=10)
    ax4.set_ylabel(r'$\alpha$', fontsize=15)
    ax4.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
    return
def plotter_noT():
    os.chdir(path)
    fig =plt.figure(figsize=[12,10],tight_layout=True)
    ax = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)

    msd = np.loadtxt('msd')
    msdx = np.loadtxt('msd_x')
    msdy = np.loadtxt('msd_y')
    ax.plot(np.mean(msd,axis=0)[n:])
    ax2.plot(np.mean(msdx,axis=0)[n:])
    ax3.plot(np.mean(msdy,axis=0)[n:])
    ax.set_xlabel('Lag time', fontsize=15)
    ax2.set_xlabel('Lag time', fontsize=15)
    ax3.set_xlabel('Lag time', fontsize=15)
    ax.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$', fontsize=15)
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$', fontsize=15)
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$', fontsize=15)
    plt.show()
    return
def plotter_2batchdata_normal():
    ########### fitted lines ###########
    t_t = np.linspace(1e-2,100,100)
    param = [0.867, 0.871, 0.869, 0.866, 0.860]  # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
    xparam = [0.867, 0.881, 0.891, 0.899, 0.908]
    yparam = [0.867, 0.858, 0.840, 0.820, 0.792]
    alphax_sem_sim = [3.24e-3, 2.45e-3, 2.05e-3, 1.79e-3,
                      1.66e-3]  # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
    alphay_sem_sim = [3.27e-3, 2.39e-3, 2.08e-3, 2.03e-3, 2.48e-3]
    alpha_sem_sim = [2.34e-3, 1.73e-3, 1.41e-3, 1.23e-3, 1.15e-3]
    # print('Mean exponent:' + str(np.mean(param[:, 1])) + '+-' + str(np.std(param[:, 1])))
    # print('X Mean exponent:' + str(np.mean(xparam[:, 1])) + '+-' + str(np.std(xparam[:, 1])))
    # print('Y Mean exponent:' + str(np.mean(yparam[:, 1])) + '+-' + str(np.std(yparam[:, 1])))
    y1 = 2e-2*(t_t/t_t[1])**0.85
    y2 = 2e-2*(t_t/t_t[1])**0.92
    y31 = 1e-2*(t_t/t_t[1])**0.95
    y32 = 4e-3*(t_t)**0.57


    fig =plt.figure(figsize=[12,10],tight_layout=True)
    ax = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    t = np.arange(1,n)
    ind = 1
    cmap = cm.get_cmap("inferno").colors
    colorint = 30
    colorshift = 0.8
    xlimr = 3000
    marker = ['o','v','^','<','>','s','+','x']
    for folder in dirlist:
        msd, msdx, msdy = msd_individualclips(path+folder)
        ax.plot(t, np.mean(msd,axis=0)[n:], color = cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=dirlist[ind-1])
        ax2.plot(t, np.mean(msdx,axis=0)[n:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=dirlist[ind-1])
        ax3.plot(t, np.mean(msdy,axis=0)[n:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=dirlist[ind-1])
        ind+=1
    # ax.plot(t_t,y1,'--',color='r', label=r'$\alpha = 0.85$')
    # ax2.plot(t_t,y2,'--',color='r', label=r'$\alpha = 0.92$')
    # ax3.plot(t_t,y31,'--',color='r', label=r'$\alpha = 0.95$')
    # ax3.plot(t_t,y32,ls='dashdot',color='b', label=r'$\alpha = 0.57$')

    # t = np.arange(0,500)
    # msd_scale = 2e-2*(t)**(0.8)
    #
    # ax.errorbar(x=np.arange(1,n), y=msd_mean[int(n):], yerr=msd_err[int(n):])
    # ax.plot(t,msd_scale,'--')
    ax.set_xlim([0, xlimr])
    ax.set_ylim([0, 1])
    ax.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$',fontsize=15)
    ax.set_xlabel('Lag time',fontsize=15)
    ax.legend(fontsize=10)
    #
    #
    # ax2.errorbar(x=np.arange(1,n), y=msdx_mean[int(n):], yerr=msdx_err[int(n):])
    # ax2.plot(t,msd_scale,'--')
    ax2.set_xlim([0, 1500])
    ax2.set_ylim([0, 2])
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$',fontsize=15)
    ax2.set_xlabel('Lag time',fontsize=15)
    ax2.legend(fontsize=10)

    #
    # ax3.errorbar(x=np.arange(1,n), y=msdy_mean[int(n):], yerr=msdy_err[int(n):])
    # ax3.plot(t,msd_scale,'--')
    ax3.set_xlim([0, 800])
    ax3.set_ylim([0, 0.5])
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$',fontsize=15)
    ax3.set_xlabel('Lag time',fontsize=15)
    ax3.legend(fontsize=10)

    x = np.arange(0,len(param))
    ax4.errorbar(x, param, yerr = alpha_sem_sim, marker='s', label=r'$\langle(\Delta r)^2\rangle$', markersize=10,capsize=4, ls='None')
    ax4.errorbar(x, xparam,yerr = alphax_sem_sim, marker='o', label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=10, capsize=4, ls='None')
    ax4.errorbar(x, yparam,yerr = alphay_sem_sim, marker='^', label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=10, capsize=4, ls='None')
    ax4.set_xticks(x)
    ax4.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax4.legend(loc = 'lower left', fontsize=10)
    ax4.set_ylabel(r'$\alpha$', fontsize=15)
    ax4.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
    return
def plotter_2batchdata_calibrated():
    ########### fitted lines ###########
    t_t = np.linspace(1e-2,100,100)
    # param = [0.867, 0.871, 0.869, 0.866, 0.860]  # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
    # xparam = [0.867, 0.881, 0.891, 0.899, 0.908] # Uncalibrated fitting for 60frames
    # yparam = [0.867, 0.858, 0.840, 0.820, 0.792]
    #
    #
    # alphax_sem_sim = [3.24e-3, 2.45e-3, 2.05e-3, 1.79e-3,
    #                   1.66e-3]  # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
    # alphay_sem_sim = [3.27e-3, 2.39e-3, 2.08e-3, 2.03e-3, 2.48e-3]
    # alpha_sem_sim = [2.34e-3, 1.73e-3, 1.41e-3, 1.23e-3, 1.15e-3]

    param = [0.69, 0.73, 0.73, 0.72, 0.73]  # Time scale calibrated fitting
    xparam = [0.68, 0.75, 0.78, 0.78, 0.81]
    yparam = [0.69, 0.71, 0.67, 0.61, 0.60]
    alpha_sem_sim = [6e-3, 3e-3, 3e-3, 3e-3, 2e-3]
    alphax_sem_sim = [6e-3, 5e-3, 4e-3, 4e-3,
                      3e-3]  # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
    alphay_sem_sim = [6e-3, 4e-3, 4e-3, 4e-3, 4e-3]
    # print('Mean exponent:' + str(np.mean(param[:, 1])) + '+-' + str(np.std(param[:, 1])))
    # print('X Mean exponent:' + str(np.mean(xparam[:, 1])) + '+-' + str(np.std(xparam[:, 1])))
    # print('Y Mean exponent:' + str(np.mean(yparam[:, 1])) + '+-' + str(np.std(yparam[:, 1])))
    y1 = 2e-2*(t_t/t_t[1])**0.85
    y2 = 2e-2*(t_t/t_t[1])**0.92
    y31 = 1e-2*(t_t/t_t[1])**0.95
    y32 = 4e-3*(t_t)**0.57


    fig =plt.figure(figsize=[12,10],tight_layout=True)
    ax = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ind = 1
    cmap = cm.get_cmap("inferno").colors
    colorint = 30
    colorshift = 0.8
    xliml = 1e-1
    xlimr = 50
    marker = ['o','v','^','<','>','s','+','x']
    path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/msdsim_calibration/Simulation/'
    os.chdir(path)
    for folder in datalist:
        os.chdir(path+folder)
        msd = np.loadtxt('msdsim_mean.txt')
        msdx = np.loadtxt('msdxsim_mean.txt')
        msdy = np.loadtxt('msdysim_mean.txt')
        l = len(msd)
        start = int((l-1)/2)
        end = len(msd[start:])
        t = np.arange(0, len(msd))
        ax.plot(t[:end]*ratio[ind-1], msd[start:], color = cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=datalist[ind-1])
        ax2.plot(t[:end]*ratio[ind-1], msdx[start:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=datalist[ind-1])
        ax3.plot(t[:end]*ratio[ind-1], msdy[start:], color=cmap[int((ind + colorshift) * colorint)],marker=marker[ind],label=datalist[ind-1])
        # ax.plot(msd[start:], color=cmap[int((ind + colorshift) * colorint)],
        #         marker=marker[ind], label=datalist[ind - 1])
        # ax2.plot(msdx[start:], color=cmap[int((ind + colorshift) * colorint)],
        #          marker=marker[ind], label=datalist[ind - 1])
        # ax3.plot(msdy[start:], color=cmap[int((ind + colorshift) * colorint)],
        #          marker=marker[ind], label=datalist[ind - 1])
        ind+=1
    # ax.plot(t_t,y1,'--',color='r', label=r'$\alpha = 0.85$')
    # ax2.plot(t_t,y2,'--',color='r', label=r'$\alpha = 0.92$')
    # ax3.plot(t_t,y31,'--',color='r', label=r'$\alpha = 0.95$')
    # ax3.plot(t_t,y32,ls='dashdot',color='b', label=r'$\alpha = 0.57$')

    # t = np.arange(0,500)
    # msd_scale = 2e-2*(t)**(0.8)
    #
    # ax.errorbar(x=np.arange(1,n), y=msd_mean[int(n):], yerr=msd_err[int(n):])
    # ax.plot(t,msd_scale,'--')
    ax.set_xlim([xliml, xlimr])
    ax.set_ylim([3e-2, 10])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$',fontsize=15)
    ax.set_xlabel('Lag time',fontsize=15)
    ax.legend(fontsize=10)
    #
    #
    # ax2.errorbar(x=np.arange(1,n), y=msdx_mean[int(n):], yerr=msdx_err[int(n):])
    # ax2.plot(t,msd_scale,'--')
    ax2.set_xlim([xliml, xlimr])
    ax2.set_ylim([3e-2, 10])
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$',fontsize=15)
    ax2.set_xlabel('Lag time',fontsize=15)
    ax2.legend(fontsize=10)

    #
    # ax3.errorbar(x=np.arange(1,n), y=msdy_mean[int(n):], yerr=msdy_err[int(n):])
    # ax3.plot(t,msd_scale,'--')
    ax3.set_xlim([xliml, xlimr])
    ax3.set_ylim([1e-2, 10])
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$',fontsize=15)
    ax3.set_xlabel('Lag time',fontsize=15)
    ax3.legend(fontsize=10)

    x = np.arange(0,len(param))
    ax4.errorbar(x, param, yerr = alpha_sem_sim, marker='s', label=r'$\langle(\Delta r)^2\rangle$', markersize=10,capsize=4, ls='None')
    ax4.errorbar(x, xparam,yerr = alphax_sem_sim, marker='o', label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=10, capsize=4, ls='None')
    ax4.errorbar(x, yparam,yerr = alphay_sem_sim, marker='^', label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=10, capsize=4, ls='None')
    ax4.set_xticks(x)
    ax4.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax4.legend(loc = 'lower left', fontsize=10)
    ax4.set_ylabel(r'$\alpha$', fontsize=15)
    ax4.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
    return
def timerescale():
    datalist = ['ecc0', 'ecc06', 'ecc08', 'ecc09', 'ecc095']
    path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/msdsim_calibration/'
    param_sim = [0.69, 0.73, 0.73, 0.72, 0.73]  # Time scale calibrated fitting
    aparam_sim = [0.022,0.019,0.018,0.019,0.018]
    os.chdir(path)
    threshold = 0.85
    ratio = []
    def indexfinder(exp_msd):
        exp_platform = exp_msd[int(len(exp_msd) / 3):int(len(exp_msd) / 3 * 2)]
        indarr = np.arange(0, len(exp_msd))
        a = indarr[exp_msd > (np.mean(exp_platform) * threshold)][0]
        b = np.mean(exp_platform)
        return a
    ind = 1
    fig = plt.figure(figsize=[10,12])
    for item in datalist:
        # loading experimental msd data
        if item == 'ecc0':
            exp_msd = np.loadtxt('./Experiment/'+item+'/'+item+'_msd.txt')[0,:]
            ftime = 0.021
            pixel = 0.11
        else:
            exp_msd = np.loadtxt('./Experiment/' + item + '/' + item + '_msd.txt') # in pixel unit
            ftime = 0.0588
            pixel = 0.16
        def func(x,a,b):
            return a*x**b
        # loading simulation msd data
        sim_msd = np.loadtxt('./Simulation/'+item+'/msdsim_mean.txt')
        sim_msd = sim_msd[int((len(sim_msd)-1)/2):]
        sim_fit = func(np.arange(0, len(sim_msd)), aparam_sim[ind-1], param_sim[ind-1])-0.015

        exp_t = indexfinder(exp_msd)*ftime # cut off frame for experiment
        sim_t = indexfinder(sim_msd) # cut off frame for simulation (reaching platform)

        r = exp_t/sim_t
        ratio.append(r)

        t_exp = np.arange(0, len(exp_msd)) * ftime
        t_sim = np.arange(0, len(sim_msd)) * r

        # if ind == 4:
        #     exp_platform = np.mean(exp_msd[int(len(exp_msd) / 3):int(len(exp_msd) / 3 * 2)])
        #     sim_platform = np.mean(sim_msd[int(len(sim_msd) / 3):int(len(sim_msd) / 3 * 2)])
        #     sim_msd = sim_msd*exp_platform/sim_platform*pixel**2
        ax = fig.add_subplot(3,2,ind)


        # ax.plot(t_exp, exp_msd*pixel**2,label = 'Experiment_'+item, marker='x')
        ax.plot(t_sim, sim_msd, label='Simulation_'+item, marker='+')
        # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t', t_sim)
        # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd', sim_msd)
        # ax.plot(t_sim, sim_fit, label='Fitted MSD')
        ax.set_xlabel('Lag time (s)', fontsize=15)
        ax.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$', fontsize=15)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-1, 1e2])
        ax.set_ylim([1e-2, 1e1])
        ax.legend()
        ind+=1
        print(item+' ratio:'+str(r)+'s/sim frame')
        print('Fitting frame for 1s:'+str(1/r))
    plt.show()
    return ratio
def timerescale_para():
    datalist = ['ecc0', 'ecc06', 'ecc08', 'ecc09', 'ecc095']
    path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/msdsim_calibration/'
    os.chdir(path)
    threshold = 0.85
    ratio = []
    def indexfinder(exp_msd):
        exp_platform = exp_msd[int(len(exp_msd) / 3):int(len(exp_msd) / 3 * 2)]
        indarr = np.arange(0, len(exp_msd))
        a = indarr[exp_msd > (np.mean(exp_platform) * threshold)][0]
        return a
    ind = 1
    fig = plt.figure(figsize=[10,12])
    for item in datalist:
        # loading experimental msd data
        if item == 'ecc0':
            exp_msd = np.loadtxt('./Experiment/'+item+'/'+item+'_msdx.txt')[0,:]
            ftime = 0.021
            pixel = 0.11
        else:
            exp_msd = np.loadtxt('./Experiment/' + item + '/' + item + '_msdx.txt') # in pixel unit
            ftime = 0.0588
            pixel = 0.16
        # loading simulation msd data
        sim_msd = np.loadtxt('./Simulation/'+item+'/msdxsim_mean.txt')
        sim_msd = sim_msd[int((len(sim_msd)-1)/2):]

        exp_t = indexfinder(exp_msd)*ftime # cut off frame for experiment
        sim_t = indexfinder(sim_msd) # cut off frame for simulation (reaching platform)

        r = exp_t/sim_t
        ratio.append(r)

        t_exp = np.arange(0, len(exp_msd)) * ftime
        t_sim = np.arange(0, len(sim_msd)) * r

        if ind == 4:
            exp_platform = np.mean(exp_msd[int(len(exp_msd) / 3):int(len(exp_msd) / 3 * 2)])
            sim_platform = np.mean(sim_msd[int(len(sim_msd) / 3):int(len(sim_msd) / 3 * 2)])
            sim_msd = sim_msd*exp_platform/sim_platform*pixel**2
        ax = fig.add_subplot(3,2,ind)
        ax.plot(t_exp, exp_msd*pixel**2,label = 'Experiment_'+item, marker='x')
        ax.plot(t_sim, sim_msd, label='Simulation_'+item, marker='+')

        ax.set_xlabel('Lag time (s)', fontsize=15)
        ax.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$', fontsize=15)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-1, 1e2])
        ax.set_ylim([1e-2, 1e1])
        ax.legend()
        ind+=1
        print(item+' ratio:'+str(r)+'s/sim frame')
        print('Fitting frame for 1s:'+str(1/r))
    plt.show()
    return ratio
# path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/msdsim_calibration/' # For calibrated MSD plot
# os.chdir(path)
path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/' # For uncalibrated MSD plot
os.chdir(path)
ratio = timerescale()
# timerescale_para()

datalist = ['ecc0', 'ecc06', 'ecc08', 'ecc09', 'ecc095']
# plotter_2batchdata()
# exponentfit_individualclips()
# plotter_2batchdata_calibrated()

# path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/ecc0'
# os.chdir(path)
# msd,msdx,msdy = msd_individualclips_ecc0newsim(path)
# np.savetxt('msd_sim', msd)
# np.savetxt('msdx_sim', msdx)
# np.savetxt('msdy_sim', msdy)