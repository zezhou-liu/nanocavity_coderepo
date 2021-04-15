import numpy as np
import matplotlib.pyplot as plt
import module
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import Normalize
import json
import os
from scipy.optimize import curve_fit
import matplotlib.colors as mcolors
from matplotlib.transforms import Affine2D
import matplotlib.patches as patches


main_path = "/media/zezhou/Seagate Expansion Drive/McGillResearch/2019Manuscript_Analysis/Analysis/subdata/datafterlinearshift/tplasmid"
cleanmode = 1 # 1-cleaned data. The roi.json and tot_file_clean.json files have been saved in ./data folder.
              # 0-raw data. Experimental data before shift to zero.
# Read files
handle1, tot_file = module.bashload(main_path)
handle1, tot_vector = module.bashvector(handle1)
handle1, tot_vec_overlay = module.bashoverlay(handle1)

# Data clean
if cleanmode == 0:
    handle1, roi = module.bashroi(handle1) # ROI selection
    handle1, tot_file_clean = module.bashclean(handle1) # Delete points ouside ROI and attach mask to handle1.
    handle1, tot_file_shift = module.bashshift(handle1) # Shift data to zero according to YOYO-3 channel
elif cleanmode == 1:
    os.chdir(main_path+'/data')
    tot_file_clean = json.load(open('tot_file_clean.json')) # data is saved in list format
    for filename in tot_file_clean:
        tot_file_clean[filename] = np.array(tot_file_clean[filename]) # Changing format to array
    handle1.tot_file_shift = tot_file_clean

# Cleaned data re-calculate
handle1, tot_vector_clean = module.bashvector(handle1, mode='clean')
handle1, tot_vec_overlay_clean = module.bashoverlay(handle1, mode='clean')
handle1, tot_pos_overlay_shift = module.bashoverlay(handle1, mode='clean', set='position')

##########################
datalist = ['ecc0', 'ecc06', 'ecc08', 'ecc09', 'ecc095', 'ecc098', 'ecc0995']
###Data setup##########
MSD_ecc0 = []
MSDx_ecc0 =[]
MSDy_ecc0 = []
MSD_ecc06 = []
MSDx_ecc06 = []
MSDy_ecc06 = []
MSD_ecc08 = []
MSDx_ecc08 = []
MSDy_ecc08 = []
MSD_ecc09 = []
MSDx_ecc09 = []
MSDy_ecc09 = []
MSD_ecc095 = []
MSDx_ecc095 = []
MSDy_ecc095 = []
MSD_ecc098 = []
MSDx_ecc098 = []
MSDy_ecc098 = []
MSD_ecc0995 = []
MSDx_ecc0995 = []
MSDy_ecc0995 = []

# Calculate MSD, MSDx and MSDy for each video
for item in handle1.tot_file_shift:
    if 'ecc0_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc0.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc0.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc0.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    if 'ecc06_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc06.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc06.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc06.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    if 'ecc08_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc08.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc08.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc08.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    if 'ecc09_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc09.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc09.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc09.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    if 'ecc095_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc095.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc095.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc095.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    if 'ecc098_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc098.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc098.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc098.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    if 'ecc0995_' in item and 'y1x' in item:
        xtemp = handle1.tot_file_shift[item]
        ytemp = handle1.tot_file_shift[item.replace('y1x', 'y1y')]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD_ecc0995.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx_ecc0995.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy_ecc0995.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))

##################### Take average of the clips #############
msd_avg = []
std_avg = []
msdx_avg = []
stdx_avg = []
msdy_avg = []
stdy_avg = []
time = 1 # fitting time
os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/msdsim_calibration')
for item in datalist:
    msd = eval('MSD_'+item)
    msdx = eval('MSDx_'+item)
    msdy = eval('MSDy_'+item)
    nclip = len(msd)
    tot_frames = 0
    nframes = []
    for index in range(nclip):
        nframes.append(msd[index].size)
        tot_frames += msd[index].size
    # arr: clips array for one eccentricity.
    def msd_average(inmsd, nframes_in):
        # nframes_in: list of frames for each clip
        arr = np.zeros((len(nframes_in), np.max(nframes_in))) # note msd is symmetric and we only process the positive(right side) value.
        for idx in range(len(inmsd)):
            arr[idx, int((np.max(nframes_in)-1)/2):int((np.max(nframes_in)-1)/2+(nframes_in[idx]+1)/2)] = inmsd[idx][int((nframes_in[idx]-1)/2):] # save the right side value.

        # sort array from min length to max length
        indarray = np.argsort(nframes_in)
        nframes_in=np.array(nframes_in)[indarray]
        arr = arr[indarray, :]
        # weight list
        msd_temp = np.zeros(int((np.max(nframes_in)+1)/2)) # mean value of msd
        std_temp = np.zeros(int((np.max(nframes_in) + 1) / 2)) # std of msd

        stdarr = arr
        stdarr[stdarr==0] = np.nan
        for idx in range(len(inmsd)):
            if idx == 0:
                msd_temp[0:int((nframes_in[idx]+1)/2)] =np.average(arr[:,int((np.max(nframes_in)-1)/2):int((np.max(nframes_in)-1)/2+(nframes_in[idx]+1)/2)], axis=0, weights=nframes_in)
                std_temp[0:int((nframes_in[idx]+1)/2)] = np.nanstd(stdarr[:,int((np.max(nframes_in)-1)/2):int((np.max(nframes_in)-1)/2+(nframes_in[idx]+1)/2)], axis=0)
            else:
                msd_temp[int((nframes_in[idx - 1]+1)/2):int((nframes_in[idx]+1)/2)] = np.average(arr[idx:,
                                                                                                 int((np.max(nframes_in)-1)/2+(nframes_in[idx - 1]+1)/2):
                                                                                                 int((np.max(nframes_in)-1)/2+(nframes_in[idx]+1)/2)], axis=0, weights = nframes_in[idx:])
                std_temp[int((nframes_in[idx - 1] + 1) / 2):int((nframes_in[idx] + 1) / 2)] = np.nanstd(stdarr[idx:,
                                                                                                         int((np.max(
                                                                                                             nframes_in) - 1) / 2 + (
                                                                                                                         nframes_in[
                                                                                                                             idx - 1] + 1) / 2):
                                                                                                         int((np.max(
                                                                                                             nframes_in) - 1) / 2 + (
                                                                                                                         nframes_in[
                                                                                                                             idx] + 1) / 2)],
                                                                                                         axis=0)

        return msd_temp, std_temp

    msd_t, std_t = msd_average(msd, nframes)
    msd_avg.append(msd_t)
    std_avg.append(std_t)
    np.savetxt(item+'_msd.txt', msd_t)
    np.savetxt(item + '_std.txt', std_t)
    msd_t, std_t = msd_average(msdx, nframes)
    msdx_avg.append(msd_t)
    stdx_avg.append(std_t)
    np.savetxt(item + '_msdx.txt', msd_t)
    np.savetxt(item + '_stdx.txt', std_t)
    msd_t, std_t = msd_average(msdy, nframes)
    msdy_avg.append(msd_t)
    stdy_avg.append(std_t)
    np.savetxt(item + '_msdy.txt', msd_t)
    np.savetxt(item + '_stdy.txt', std_t)


#################### exponent fitting ##################
def exponentfit():
    time = 1 # fitting time (sec)
    ftime = int(time*17) # fitting frames
    def fitfunc(x,a,b):# fitting function for MSD. Exponent
        return a*x**b
    sol = []
    cov = []
    for item in msd_avg:
        msd_tmp = item[:ftime]
        t_tmp = np.linspace(5e-2,time,ftime)
        popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp, p0=[1,0.5])
        sol.append(popt)
        cov.append(pcov)
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/MSDfit/')
    print(sol)
    print(cov)
    f = open('msdfit.txt', 'a')
    np.savetxt(f,sol)
    f.close()
    f = open('msdfit_cov.txt', 'a')
    for item in cov:
        np.savetxt(f,item)
    f.close()

def exponentfit_individualvid():
    # Diffusion exponent fitting based on each video. This way is supposed to extract the errorbar.
    ftime = int(time * 17)  # fitting frames

    def fitfunc(x, a, b):  # fitting function for MSD. Exponent
        return a * x ** b

    for item in datalist:
        msd = eval('MSD_' + item)
        msdx = eval('MSDx_' + item)
        msdy = eval('MSDy_' + item)

        sol = []
        cov = []
        solx = []
        covx = []
        soly = []
        covy = []
        for clip_ind in range(len(msd)):
            msd_tmp = msd[clip_ind]
            start = int((len(msd_tmp)-1)/2) # take the right side of MSD since it's symmetric
            end = start+ftime
            msd_tmp = msd_tmp[start:end]

            msdx_tmp = msdx[clip_ind]
            start = int((len(msdx_tmp) - 1) / 2)  # take the right side of MSD since it's symmetric
            end = start + ftime
            msdx_tmp = msdx_tmp[start:end]

            msdy_tmp = msdy[clip_ind]
            start = int((len(msdy_tmp) - 1) / 2)  # take the right side of MSD since it's symmetric
            end = start + ftime
            msdy_tmp = msdy_tmp[start:end]

            t_tmp = np.linspace(0, time, ftime)
            popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp, p0=[1, 0.5])
            poptx, pcovx = curve_fit(fitfunc, t_tmp, msdx_tmp, p0=[1, 0.5])
            popty, pcovy = curve_fit(fitfunc, t_tmp, msdy_tmp, p0=[1, 0.5])

            sol.append(popt)
            cov.append(pcov)
            solx.append(poptx)
            covx.append(pcovx)
            soly.append(popty)
            covy.append(pcovy)
        os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/MSDfit/test')
        f = open('msdfit_clip_'+item+'.txt', 'a')
        np.savetxt(f, sol)
        f.close()
        f = open('msdfit_cov_clip_'+item+'.txt', 'a')
        for item11 in cov:
            np.savetxt(f, item11)
        f.close()

        f = open('msdfitx_clip_' + item + '.txt', 'a')
        np.savetxt(f, solx)
        f.close()
        f = open('msdfitx_cov_clip_' + item + '.txt', 'a')
        for item11 in covx:
            np.savetxt(f, item11)
        f.close()

        f = open('msdfity_clip_' + item + '.txt', 'a')
        np.savetxt(f, soly)
        f.close()
        f = open('msdfity_cov_clip_' + item + '.txt', 'a')
        for item11 in covy:
            np.savetxt(f, item11)
        f.close()
    return

# exponentfit_individualvid()
# load fitted value
os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/MSDfit')
# param = np.loadtxt('msdfit.txt')
cov = np.loadtxt('msdfit_cov.txt')
# xparam = np.loadtxt('msdfitx.txt')
xcov = np.loadtxt('msdfitx_cov.txt')
# yparam = np.loadtxt('msdfity.txt')
ycov = np.loadtxt('msdfity_cov.txt')

# alphax_sem = [0.020,0.020,0.038,0.036,0.015,0.024,0.034] # check msd_alphaerror.py and the summary in /MSDfit/
# alphay_sem = [0.015,0.018,0.017,0.020,0.021,0.017,0.018] # the starting time was set to 5e-2 for unknown reason. I forget why I did that
# alpha_sem = [0.021,0.018,0.028,0.031,0.013,0.023,0.030]

param = np.array([0.64,0.64,0.64,0.58,0.66,0.73,0.61])
alpha_sem = np.array([0.02,0.02,0.03,0.07,0.01,0.02,0.05])
xparam = np.array([0.62,0.67,0.7,0.63,0.75,0.82,0.69])
yparam = np.array([0.67,0.60,0.55,0.46,0.35,0.35,0.22])
alphax_sem = np.array([0.02,0.02,0.04,0.08,0.02,0.02,0.03]) # check msd_alphaerror.py and the summary in /MSDfit/. start fitting time was set to 0.
alphay_sem = np.array([0.02,0.02,0.02,0.06,0.03,0.02,0.02])

# alphax_sem_sim = [3.24e-3, 2.45e-3, 2.05e-3,1.79e-3, 1.66e-3] # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
# alphay_sem_sim = [3.27e-3, 2.39e-3, 2.08e-3, 2.03e-3, 2.48e-3]
# param_sim = [0.867,0.871,0.869,0.866,0.860]
# xparam_sim = [0.867,0.881,0.891,0.899,0.908]
# yparam_sim = [0.867,0.858,0.840, 0.820, 0.792]
#
# param_sim = [0.79,0.79,0.78,0.77,0.77] # hotfix. 1s~60 frams
# xparam_sim = [0.79,0.81,0.82,0.83,0.84]
# yparam_sim = [0.79,0.76,0.72, 0.69, 0.63]

param_sim = [0.69,0.73,0.73,0.72,0.73] # Time scale calibrated fitting
xparam_sim = [0.68,0.75,0.78,0.78,0.81]
yparam_sim = [0.69,0.71,0.67,0.61,0.60]
alphax_sem_sim = [6e-3,5e-3,4e-3,4e-3,3e-3] # check msd_alphaerror_sim.py and the summary in /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/
alphay_sem_sim = [6e-3,4e-3,4e-3,4e-3,4e-3]

os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/ecc0/') # The path for the additional data

msd_add_ecc0 = np.loadtxt('msd_mean') # additional msd
msdx_add_ecc0 = np.loadtxt('msdx_mean')
msdy_add_ecc0 = np.loadtxt('msdy_mean')

idx = 0
err = []
xerr = []
yerr = []
fontsize = 18

colorint = 30
colorshift = 0.8
ecc0color = int((1 + colorshift) * colorint)
t_additional = np.arange(0,len(msd_add_ecc0[0,:]))*0.021
# for idx in range(np.shape(cov)[0]):
#     if idx%2!=0:
#         err.append(cov[idx,1]**2)
#         xerr.append(xcov[idx, 1] ** 2)
#         yerr.append(ycov[idx, 1] ** 2)
#
# print('Mean exponent:'+str(np.mean(param[:,1]))+'+-'+str(np.std(param[:,1])))
# print('X Mean exponent:'+str(np.mean(xparam[:,1]))+'+-'+str(np.std(xparam[:,1])))
# print('Y Mean exponent:'+str(np.mean(yparam[:,1]))+'+-'+str(np.std(yparam[:,1])))
#################### Log Plot ##########################3
def logplot():
    fig = plt.figure(figsize=[12,10],tight_layout=True)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    marker = ['o','v','^','<','>','s','+','x']
    cmap = cm.get_cmap("inferno").colors
    def plotter(ax, inmsd_avg):
        ind = 1
        cmap = cm.get_cmap("inferno").colors
        for item in inmsd_avg:
            if ind == 1:
                ind += 1
                continue
            x = np.arange(0,len(item))/17
            ax.plot(x, item/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=datalist[ind-1])
            ind +=1
        ax.set_xlabel('Lag time (s)', fontsize=15)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
        ax.set_xlim([5e-2,300])
        ax.set_ylim([3e-2,10])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()

    x = np.linspace(5e-2,50,2000)
    y = 1e-1*(x/(5e-2))**(0.82)
    y1 = (8e-2)*(x/(5e-2))**(0.78)
    y3 = (5e-2)*(x/(5e-2))**(0.84)
    y4 = (1e-2)*(x/(5e-2))**(0.34)
    ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    ax2.plot(x,y,'--', label=r'$\alpha = 0.82$',color='r')
    ax3.plot(x,y3,'--', label=r'$\alpha = 0.84$',color='r')
    ax3.plot(x, y4, '-.', label=r'$\alpha = 0.34$', color='b')
    ax1.plot(t_additional, msd_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
             label='ecc0')
    plotter(ax1, msd_avg)

    ax1.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$', fontsize=15)
    ax2.plot(t_additional, msdx_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
             label='ecc0')
    plotter(ax2, msdx_avg)

    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$', fontsize=15)
    ax3.plot(t_additional, msdy_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
             label='ecc0')
    plotter(ax3,msdy_avg)

    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$', fontsize=15)
    ax3.set_xlim([6e-2, 300])
    ax3.set_ylim([1e-2, 10])
    ax3.set_xlim([6e-2, 100])
    ax2.set_xlim([6e-2, 100])
    ax1.set_xlim([6e-2, 100])
    ax1.legend(loc='upper right')

    ax = fig.add_subplot(224)
    x = np.arange(0,len(param))
    ax.errorbar(x, param, yerr=np.array(alphax_sem), marker='^', label=r'$\langle(\Delta r)^2\rangle$', markersize=10, capsize=5, ls='None')
    ax.errorbar(x, xparam, yerr=np.array(alphax_sem), marker='<',
                label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=8, capsize=5, ls='None')
    ax.errorbar(x, yparam, yerr=np.array(alphay_sem), marker='>',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    ax.set_xticks(x)
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95', '0.98', '0.995'])
    ax.legend(loc = 'lower left', fontsize=15)
    ax.set_ylabel(r'$\alpha$', fontsize=15)
    ax.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
def logplot_v2():
    fig = plt.figure(figsize=[6.5,6.5*10/6])
    # ax1 = fig.add_subplot(221)

    ax3 = fig.add_subplot(211)
    ax2 = fig.add_axes([0.55, 0.7, 0.35, 0.17])
    marker = ['o','v','^','<','>','s','+','x']
    def plotter(ax, inmsd_avg):
        ind = 1
        cmap = cm.get_cmap("inferno").colors
        colorint = 30
        colorshift = 0.8
        for item in inmsd_avg:
            if 0 == 0:
                x = np.arange(0,len(item))/17
                ax.plot(x, item/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=datalist[ind-1])
            ind +=1
        ax.set_xlabel('Lag time (s)', fontsize=15)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
        ax.set_xlim([5e-2,300])
        ax.set_ylim([3e-2,10])
        ax.set_xscale('log')
        ax.set_yscale('log')

    x = np.linspace(5e-2,50,2000)
    y = 1e-1*(x/(5e-2))**(0.82)
    y1 = (8e-2)*(x/(5e-2))**(0.78)
    y3 = (5e-2)*(x/(5e-2))**(0.84)
    y4 = (1e-2)*(x/(5e-2))**(0.34)
    # ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    ax2.plot(x,y,'--', label=r'$\alpha = 0.82$',color='r')
    ax3.plot(x,y3,'--', label=r'$\alpha = 0.84$',color='r')
    ax3.plot(x, y4, '-.', label=r'$\alpha = 0.34$', color='b')
    # plotter(ax1, msd_avg)
    # ax1.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
    plotter(ax2, msdx_avg)
    ax2.set_xlabel('Lag time (s)', fontsize=13)
    ax2.set_ylabel(r'MSD$_{major}$ $(\mu m^2)$', fontsize=13)
    plotter(ax3,msdy_avg)
    ax3.set_ylabel(r'MSD$_{minor}$ $(\mu m^2)$', fontsize=15)
    ax3.legend(fontsize=10)
    ax3.set_xlim([6e-2, 300])
    ax3.set_ylim([1e-2, 10])
    # ax1.legend(loc='upper right')

    ax = fig.add_subplot(212)
    x = np.arange(0,len(param))
    # ax.plot(x, param[:,1],'s', label='MSD', markersize=10)
    ax.plot(x, xparam,'o',color='k', label='MSD major', markersize=10)
    ax.plot(x, yparam,'p',color='r', label='MSD minor', markersize=10)

    param_sim = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msd/msdfit.txt')
    xparam_sim = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdx/msdfit.txt')
    yparam_sim = np.loadtxt(
        '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdy/msdfit.txt')

    ax.set_xticklabels(['', '0', '0.6', '0.8', '0.9', '0.95', '0.98', '0.995'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=15)
    ax.set_xlabel('Eccentricity', fontsize=15)

    ax =  fig.add_axes([0.24, 0.13, 0.3, 0.13])
    # ax.plot(x, param_sim[:, 1], 's', label='MSD simulation', markersize=10)
    ax.plot(x, xparam_sim, 'o',color='k', label='MSD major simulation', markersize=10)
    ax.plot(x, yparam_sim, 'p',color='r', label='MSD minor simulation', markersize=10)
    ax.set_xticks([0,2,4,6])
    ax.set_xticklabels(['0', '0.8',  '0.95', '0.995'])
    # ax.legend(loc = 'lower left', fontsize=10)
    ax.set_ylabel(r'$\alpha^{MC}$', fontsize=15)
    # ax.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
def logplot_v3():
    fig = plt.figure(figsize=[8,8], tight_layout=True)
    # ax1 = fig.add_subplot(221)

    ax3 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    marker = ['o','v','^','<','>','s','+','x']
    def plotter(ax, inmsd_avg,label=0):
        ind = 1
        cmap = cm.get_cmap("inferno").colors
        colorint = 30
        colorshift = 0.8
        for item in inmsd_avg:
            if True:
                x = np.arange(0,len(item))/17
                if label==0:
                    ax.plot(x, item/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind])
                else:
                    ax.plot(x, item/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind], label = datalist[ind - 1])
            ind +=1
        ax.set_xlabel('Lag time (s)', fontsize=fontsize)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=fontsize)
        ax.set_xlim([5e-2,300])
        ax.set_ylim([3e-2,10])
        ax.set_xscale('log')
        ax.set_yscale('log')

    x = np.linspace(5e-2,50,2000)
    y = 1e-1*(x/(5e-2))**(0.82)
    y1 = (8e-2)*(x/(5e-2))**(0.78)
    y3 = (5e-2)*(x/(5e-2))**(0.84)
    y4 = (1e-2)*(x/(5e-2))**(0.34)
    # ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    ax2.plot(x,y,'--', label=r'$\alpha = 0.82$',color='r')
    ax3.plot(x,y3,'--', label=r'$\alpha = 0.84$',color='r')
    ax3.plot(x, y4, '-.', label=r'$\alpha = 0.34$', color='b')
    # plotter(ax1, msd_avg)
    # ax1.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
    plotter(ax2, msdx_avg, label=1)
    ax2.set_xlabel('Lag time (s)', fontsize=fontsize)
    ax2.set_ylabel(r'MSD$_{major}$ $(\mu m^2)$', fontsize=fontsize)
    plotter(ax3,msdy_avg)
    ax2.legend(fontsize=9)
    ax3.set_ylabel(r'MSD$_{minor}$ $(\mu m^2)$', fontsize=fontsize)
    ax3.legend(fontsize=9)
    ax3.set_xlim([6e-2, 300])
    ax3.set_ylim([1e-2, 10])
    # ax1.legend(loc='upper right')

    ax = fig.add_subplot(223)
    x = np.arange(0,len(param[:,1]))
    # ax.plot(x, param[:,1],'s', label='MSD', markersize=10)
    ax.plot(x, xparam[:,1],'o',color='k', label='MSD major', markersize=10)
    ax.plot(x, yparam[:,1],'p',color='r', label='MSD minor', markersize=10)

    param_sim = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msd/msdfit.txt')
    xparam_sim = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdx/msdfit.txt')
    yparam_sim = np.loadtxt(
        '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdy/msdfit.txt')
    ax.set_xticks(x)
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95', '0.98', '0.995'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)

    ax =  fig.add_subplot(224)
    # ax.plot(x, param_sim[:, 1], 's', label='MSD simulation', markersize=10)
    ax.plot(x[:-1], xparam_sim[:-1, 1], 'o',color='k', label='MSD major simulation', markersize=10)
    ax.plot(x[:-1], yparam_sim[:-1, 1], 'p',color='r', label='MSD minor simulation', markersize=10)
    ax.set_xticks(x[:-1])
    ax.set_xticklabels(['0', '0.6','0.8', '0.9','0.95', '0.98'])
    # ax.legend(loc = 'lower left', fontsize=10)
    ax.set_ylabel(r'$\alpha^{MC}$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    plt.show()
def logplot_v4():
    fig = plt.figure(figsize=[8,8], tight_layout = True)
    # ax1 = fig.add_subplot(221)

    ax3 = fig.add_subplot(222)
    ax2 = fig.add_subplot(221)
    marker = ['o','v','^','<','>','s','+','x']
    def plotter(ax, inmsd_avg, label=0):
        ind = 1
        cmap = cm.get_cmap("inferno").colors
        colorint = 30
        colorshift = 0.8
        for indxx in range(len(inmsd_avg)):
            if True:
                x = np.arange(0,len(inmsd_avg[indxx]))/17
                if label==0:
                    ax.plot(x, inmsd_avg[indxx]/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind])
                else:
                    ax.plot(x, inmsd_avg[indxx]/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind], label = datalist[ind - 1])
            ind +=1
        ax.set_xlabel('Lag time (s)', fontsize=fontsize)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=fontsize)
        ax.set_xlim([5e-2,300])
        ax.set_ylim([3e-2,10])
        ax.set_xscale('log')
        ax.set_yscale('log')

    x = np.linspace(5e-2,50,2000)
    y = 1e-1*(x/(5e-2))**(0.82)
    y1 = (8e-2)*(x/(5e-2))**(0.78)
    y3 = (5e-2)*(x/(5e-2))**(0.84)
    y4 = (1e-2)*(x/(5e-2))**(0.34)
    # ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    ax2.plot(x,y,'--', label=r'$\alpha = 0.82$',color='r')
    ax3.plot(x,y3,'--', label=r'$\alpha = 0.84$',color='r')
    ax3.plot(x, y4, '-.', label=r'$\alpha = 0.34$', color='b')
    # plotter(ax1, msd_avg)
    # ax1.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
    plotter(ax2, msdx_avg, label=1)
    ax2.set_xlabel('Lag time (s)', fontsize=fontsize)
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$', fontsize=fontsize)
    plotter(ax3,msdy_avg)
    ax2.legend(fontsize=9)
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$', fontsize=fontsize)
    ax3.legend(fontsize=9)
    ax3.set_xlim([6e-2, 100])
    ax2.set_xlim([6e-2, 100])
    ax3.set_ylim([1e-2, 10])
    # ax2.set_ylim([1e-2, 10])
    # ax1.legend(loc='upper right')

    ax = fig.add_subplot(223)
    x = np.arange(0,len(param))
    # ax.plot(x, param[:,1],'s', label='MSD', markersize=10)
    # a = x[:-2]
    # b = xparam[:-2,1]
    # c = alphax_std[:-2]
    ax.errorbar(x[:-2], xparam[:-2], yerr=np.array(alphax_sem[:-2]),marker='<',color='r', label= r'$\langle(\Delta r)^2\rangle_{||}$', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], yparam[:-2,1], yerr=np.array(alphay_sem[:-2]),marker='>',color='r', label= r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], xparam_sim, yerr = alphax_sem_sim, marker ='^',color='k', label=r'$\langle(\Delta r)^2\rangle_{||}$ sim', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], yparam_sim, yerr = alphay_sem_sim, marker ='v',color='k', label=r'$\langle(\Delta r)^2\rangle_{\perp}$ sim', markersize=8, capsize=5,ls='None')

    ax.set_xticks(x[:-2])
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    ax.set_ylim([0.5,1])

    ax = fig.add_subplot(224)
    x = np.arange(0, len(param[:, 1]))
    # ax.plot(x, param[:,1],'s', label='MSD', markersize=10)
    # a = x[:-2]
    # b = xparam[:-2,1]
    # c = alphax_std[:-2]
    # ax.errorbar(x[:-2], xparam[:-2, 1], yerr=np.array(alphax_sem[:-2]), marker='<', color='r',
                # label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], yparam[:-2], yerr=np.array(alphay_sem[:-2]), marker='>', color='r',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], xparam_sim, yerr=alphax_sem_sim, marker='^', color='k',
                # label=r'$\langle(\Delta r)^2\rangle_{||}$ sim', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], yparam_sim, yerr=alphay_sem_sim, marker='v', color='k',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$ sim', markersize=8, capsize=5, ls='None')

    ax.set_xticks(x[:-2])
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    ax.set_ylim([0.5, 1])
    # ax =  fig.add_subplot(223)
    # # ax.plot(x, param_sim[:, 1], 's', label='MSD simulation', markersize=10)
    # ax.plot(x[:-2], xparam_sim[:-2, 1], '<',color='k', label=r'$<(\Delta \mathbf{r})^2>_{||}$ sim', markersize=10)
    # ax.plot(x[:-2], yparam_sim[:-2, 1], 'p',color='r', label=r'$<(\Delta \mathbf{r})^2>_{\perp}$ sim', markersize=10)
    # ax.set_xticks(x[:-2])
    # ax.set_xticklabels(['0', '0.6','0.8', '0.9','0.95'])
    # # ax.legend(loc = 'lower left', fontsize=10)
    # ax.set_ylabel(r'$\alpha^{MC}$', fontsize=fontsize)
    # ax.set_xlabel('Eccentricity', fontsize=fontsize)
    plt.show()
def logplot_v5():
    fig = plt.figure(figsize=[8, 8], tight_layout=True)
    # ax1 = fig.add_subplot(221)

    ax3 = fig.add_subplot(222) # ydirection
    ax2 = fig.add_subplot(221) # xdirection
    marker = ['o', 'v', '^', '<', '>', 's', '+', 'x']
    # ecc0color = int((1 + colorshift) * colorint)
    # t_additional = np.arange(0,len(msd_add_ecc0[0,:]))*0.021
    cmap = cm.get_cmap("inferno").colors
    ax3.plot(t_additional, msdy_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color])

    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t', t_additional)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd',
    #            msdy_add_ecc0[0, :] * 0.105 ** 2)
    # print('a')
    def plotter(ax, inmsd_avg, label=0, flag=0):
        ind = 2
        cmap = cm.get_cmap("inferno").colors
        if flag == 1:
            ax2.plot(t_additional, msdx_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color], label='ecc0')
            # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t', t_additional)
            # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd',
            #            msdx_add_ecc0[0, :] * 0.105 ** 2)
            # print('a')
        for indxx in range(len(inmsd_avg)-1):
            if True:
                x = np.arange(0, len(inmsd_avg[indxx+1])) / 17
                if label == 0:
                    ax.plot(x, inmsd_avg[indxx+1] / 6.25 ** 2, ls='-', color=cmap[int((ind + colorshift) * colorint)],
                            marker=marker[ind])
                    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t',x)
                    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd', inmsd_avg[indxx+1] / 6.25 ** 2)
                    # print('a')
                else:
                    ax.plot(x, inmsd_avg[indxx+1] / 6.25 ** 2, ls='-', color=cmap[int((ind + colorshift) * colorint)],
                            marker=marker[ind], label=datalist[ind - 1])
                    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t',x)
                    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd', inmsd_avg[indxx+1] / 6.25 ** 2)
                    # print('a')
            ind += 1
        ax.set_xlabel('Lag time (s)', fontsize=fontsize)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=fontsize)
        ax.set_xlim([5e-2, 300])
        ax.set_ylim([3e-2, 10])
        ax.set_xscale('log')
        ax.set_yscale('log')

    x = np.linspace(5e-2, 50, 2000)
    y = 1e-1 * (x / (5e-2)) ** (0.7)
    y1 = (8e-2) * (x / (5e-2)) ** (0.78)
    y3 = (1e-1) * (x / (5e-2)) ** (0.7)
    y4 = (1e-2) * (x / (5e-2)) ** (0.3)
    # ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    ax2.plot(x, y, '--', label=r'$\alpha = 0.7$', color='r')
    ax3.plot(x, y3, '--', label=r'$\alpha = 0.7$', color='r')
    ax3.plot(x, y4, '-.', label=r'$\alpha = 0.3$', color='b')
    # plotter(ax1, msd_avg)
    # ax1.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
    plotter(ax2, msdx_avg, label=1, flag=1)

    ax2.set_xlabel('Lag time (s)', fontsize=fontsize)
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$', fontsize=fontsize)
    plotter(ax3, msdy_avg)
    ax2.legend(fontsize=9)
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$', fontsize=fontsize)
    ax3.legend(fontsize=9)
    ax3.set_xlim([6e-2, 100])
    ax2.set_xlim([6e-2, 100])
    ax3.set_ylim([1e-2, 10])
    # ax2.set_ylim([1e-2, 10])
    # ax1.legend(loc='upper right')

    ax = fig.add_subplot(223)
    x = np.arange(0, len(param))
    # ax.plot(x, param[:,1],'s', label='MSD', markersize=10)
    # a = x[:-2]
    # b = xparam[:-2,1]
    # c = alphax_std[:-2]
    ax.errorbar(x[:-2], xparam[:-2], yerr=np.array(alphax_sem[:-2]), marker='<', color='r',
                label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], yparam[:-2,1], yerr=np.array(alphay_sem[:-2]),marker='>',color='r', label= r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], xparam_sim, yerr=alphax_sem_sim, marker='^', color='k',
                label=r'$\langle(\Delta r)^2\rangle_{||}$ sim', markersize=8, capsize=5, ls='None')

    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/para', xparam[:-2])
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/paraerr', np.array(alphax_sem[:-2]))
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simpara', xparam_sim)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simparaerr',
    #            alphax_sem_sim)
    # ax.errorbar(x[:-2], yparam_sim, yerr = alphay_sem_sim, marker ='v',color='k', label=r'$\langle(\Delta r)^2\rangle_{\perp}$ sim', markersize=8, capsize=5,ls='None')

    ax.set_xticks(x[:-2])
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    ax.set_ylim([0.33, 1])

    ax = fig.add_subplot(224)
    x = np.arange(0, len(param))

    ax.errorbar(x[:-2], yparam[:-2], yerr=np.array(alphay_sem[:-2]), marker='>', color='r',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')

    # ax.errorbar(x[:-2], xparam_sim, yerr=alphax_sem_sim, marker='^', color='k',
    # label=r'$\langle(\Delta r)^2\rangle_{||}$ sim', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], yparam_sim, yerr=alphay_sem_sim, marker='v', color='k',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$ sim', markersize=8, capsize=5, ls='None')
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/para', yparam[:-2])
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/paraerr',
    #            np.array(alphay_sem[:-2]))
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simpara', yparam_sim)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simparaerr',
    #            alphay_sem_sim)
    ax.set_xticks(x[:-2])
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    ax.set_ylim([0.33, 1])
    # ax =  fig.add_subplot(223)
    # # ax.plot(x, param_sim[:, 1], 's', label='MSD simulation', markersize=10)
    # ax.plot(x[:-2], xparam_sim[:-2, 1], '<',color='k', label=r'$<(\Delta \mathbf{r})^2>_{||}$ sim', markersize=10)
    # ax.plot(x[:-2], yparam_sim[:-2, 1], 'p',color='r', label=r'$<(\Delta \mathbf{r})^2>_{\perp}$ sim', markersize=10)
    # ax.set_xticks(x[:-2])
    # ax.set_xticklabels(['0', '0.6','0.8', '0.9','0.95'])
    # # ax.legend(loc = 'lower left', fontsize=10)
    # ax.set_ylabel(r'$\alpha^{MC}$', fontsize=fontsize)
    # ax.set_xlabel('Eccentricity', fontsize=fontsize)
    plt.show()
def logplot_supp():
    fig = plt.figure(figsize=[12,10],tight_layout=True)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    marker = ['o','v','^','<','>','s','+','x']
    cmap = cm.get_cmap("inferno").colors
    def plotter(ax, inmsd_avg):
        ind = 1
        cmap = cm.get_cmap("inferno").colors
        for item in inmsd_avg:
            if ind == 1:
                ind += 1
                continue
            x = np.arange(0,len(item))/17
            ax.plot(x, item/6.25**2, ls='-', color=cmap[int((ind+colorshift)*colorint)],marker=marker[ind],label=datalist[ind-1])
            # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t', x)
            # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd',
            #            item/6.25**2)
            # print('a')
            ind +=1
        ax.set_xlabel('Lag time (s)', fontsize=15)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
        ax.set_xlim([5e-2,300])
        ax.set_ylim([3e-2,10])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()

    x = np.linspace(5e-2,50,2000)
    y = 1e-1*(x/(5e-2))**(0.82)
    y1 = (8e-2)*(x/(5e-2))**(0.78)
    y3 = (5e-2)*(x/(5e-2))**(0.84)
    y4 = (1e-2)*(x/(5e-2))**(0.34)
    ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    ax2.plot(x,y,'--', label=r'$\alpha = 0.82$',color='r')
    ax3.plot(x,y3,'--', label=r'$\alpha = 0.84$',color='r')
    ax3.plot(x, y4, '-.', label=r'$\alpha = 0.34$', color='b')
    ax1.plot(t_additional, msd_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
             label='ecc0')
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/t',t_additional)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/msd', msd_add_ecc0[0, :] * 0.105 ** 2)
    # print('a')
    plotter(ax1, msd_avg)

    ax1.set_ylabel(r'$\langle(\Delta r)^2\rangle$ $(\mu m^2)$', fontsize=15)
    ax2.plot(t_additional, msdx_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
             label='ecc0')
    plotter(ax2, msdx_avg)

    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$', fontsize=15)
    ax3.plot(t_additional, msdy_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
             label='ecc0')
    plotter(ax3,msdy_avg)

    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$', fontsize=15)
    ax3.set_xlim([6e-2, 300])
    ax3.set_ylim([1e-2, 10])
    ax3.set_xlim([6e-2, 100])
    ax2.set_xlim([6e-2, 100])
    ax1.set_xlim([6e-2, 100])
    ax1.legend(loc='upper right')

    ax = fig.add_subplot(224)
    x = np.arange(0,len(param))
    ax.errorbar(x, param, yerr=np.array(alphax_sem), marker='^', label=r'$\langle(\Delta r)^2\rangle$', markersize=10, capsize=5, ls='None')
    ax.errorbar(x, xparam, yerr=np.array(alphax_sem), marker='<',
                label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=8, capsize=5, ls='None')
    ax.errorbar(x, yparam, yerr=np.array(alphay_sem), marker='>',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/para', param)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/parasem', np.array(alphax_sem))
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xpara', xparam)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xparasem',
    #            np.array(alphax_sem))
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/ypara', yparam)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/yparasem',
    #            np.array(alphay_sem))
    ax.set_xticks(x)
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95', '0.98', '0.995'])
    ax.legend(loc = 'lower left', fontsize=15)
    ax.set_ylabel(r'$\alpha$', fontsize=15)
    ax.set_xlabel('Eccentricity', fontsize=15)
    plt.show()
def normalplot():
    fig = plt.figure(figsize=[8, 8], tight_layout=True)
    # ax1 = fig.add_subplot(221)

    ax3 = fig.add_subplot(222)  # ydirection
    ax2 = fig.add_subplot(221)  # xdirection
    marker = ['o', 'v', '^', '<', '>', 's', '+', 'x']
    # ecc0color = int((1 + colorshift) * colorint)
    # t_additional = np.arange(0,len(msd_add_ecc0[0,:]))*0.021
    cmap = cm.get_cmap("inferno").colors
    ax3.plot(t_additional, msdy_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color])

    def plotter(ax, inmsd_avg, label=0, flag=0):
        ind = 2
        cmap = cm.get_cmap("inferno").colors
        if flag == 1:
            ax2.plot(t_additional, msdx_add_ecc0[0, :] * 0.105 ** 2, marker=marker[1], color=cmap[ecc0color],
                     label='ecc0')
        for indxx in range(len(inmsd_avg) - 3):
            if True:
                x = np.arange(0, len(inmsd_avg[indxx + 1])) / 17
                if label == 0:
                    ax.plot(x, inmsd_avg[indxx + 1] / 6.25 ** 2, ls='-', color=cmap[int((ind + colorshift) * colorint)],
                            marker=marker[ind])
                else:
                    ax.plot(x, inmsd_avg[indxx + 1] / 6.25 ** 2, ls='-', color=cmap[int((ind + colorshift) * colorint)],
                            marker=marker[ind], label=datalist[ind - 1])
            ind += 1
        ax.set_xlabel('Lag time (s)', fontsize=fontsize)
        ax.set_ylabel(r'MSD $(\mu m^2)$', fontsize=fontsize)
        ax.set_xlim([5e-2, 300])
        ax.set_ylim([3e-2, 10])
        # ax.set_xscale('log')
        # ax.set_yscale('log')

    x = np.linspace(5e-2, 50, 2000)
    # y = 1e-1 * (x / (5e-2)) ** (0.82)
    y1 = (8e-2) * (x / (5e-2)) ** (0.78)
    y3 = (5e-2) * (x / (5e-2)) ** (0.84)
    y4 = (1e-2) * (x / (5e-2)) ** (0.34)
    # # ax1.plot(x,y1,'--', label=r'$\alpha = 0.78$',color='r')
    # ax2.plot(x, y, '--', label=r'$\alpha = 0.82$', color='r')
    # ax3.plot(x, y3, '--', label=r'$\alpha = 0.84$', color='r')
    # ax3.plot(x, y4, '-.', label=r'$\alpha = 0.34$', color='b')
    # plotter(ax1, msd_avg)
    # ax1.set_ylabel(r'MSD $(\mu m^2)$', fontsize=15)
    plotter(ax2, msdx_avg, label=1, flag=1)

    ax2.set_xlabel('Lag time (s)', fontsize=fontsize)
    ax2.set_ylabel(r'$\langle(\Delta r)^2\rangle_{||}$ $(\mu m^2)$', fontsize=fontsize)
    plotter(ax3, msdy_avg)
    ax2.legend(fontsize=9)
    ax3.set_ylabel(r'$\langle(\Delta r)^2\rangle_{\perp}$ $(\mu m^2)$', fontsize=fontsize)
    ax3.legend(fontsize=9)
    ax3.set_xlim([0, 5])
    ax2.set_xlim([0, 30])
    ax2.set_ylim([0,2])
    ax3.set_ylim([0, 0.5])
    # ax2.set_ylim([1e-2, 10])
    # ax1.legend(loc='upper right')

    ax = fig.add_subplot(223)
    x = np.arange(0, len(param))
    # ax.plot(x, param[:,1],'s', label='MSD', markersize=10)
    # a = x[:-2]
    # b = xparam[:-2,1]
    # c = alphax_std[:-2]
    ax.errorbar(x[:-2], xparam[:-2], yerr=np.array(alphax_sem[:-2]), marker='<', color='r',
                label=r'$\langle(\Delta r)^2\rangle_{||}$', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], yparam[:-2,1], yerr=np.array(alphay_sem[:-2]),marker='>',color='r', label= r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], xparam_sim, yerr=alphax_sem_sim, marker='^', color='k',
                label=r'$\langle(\Delta r)^2\rangle_{||}$ sim', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], yparam_sim, yerr = alphay_sem_sim, marker ='v',color='k', label=r'$\langle(\Delta r)^2\rangle_{\perp}$ sim', markersize=8, capsize=5,ls='None')

    ax.set_xticks(x[:-2])
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    ax.set_ylim([0.33, 1])

    ax = fig.add_subplot(224)
    x = np.arange(0, len(param))

    ax.errorbar(x[:-2], yparam[:-2], yerr=np.array(alphay_sem[:-2]), marker='>', color='r',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$', markersize=8, capsize=5, ls='None')
    # ax.errorbar(x[:-2], xparam_sim, yerr=alphax_sem_sim, marker='^', color='k',
    # label=r'$\langle(\Delta r)^2\rangle_{||}$ sim', markersize=8, capsize=5, ls='None')
    ax.errorbar(x[:-2], yparam_sim, yerr=alphay_sem_sim, marker='v', color='k',
                label=r'$\langle(\Delta r)^2\rangle_{\perp}$ sim', markersize=8, capsize=5, ls='None')

    ax.set_xticks(x[:-2])
    ax.set_xticklabels(['0', '0.6', '0.8', '0.9', '0.95'])
    ax.legend()
    ax.set_ylabel(r'$\alpha$', fontsize=fontsize)
    ax.set_xlabel('Eccentricity', fontsize=fontsize)
    ax.set_ylim([0.33, 1])
    # ax =  fig.add_subplot(223)
    # # ax.plot(x, param_sim[:, 1], 's', label='MSD simulation', markersize=10)
    # ax.plot(x[:-2], xparam_sim[:-2, 1], '<',color='k', label=r'$<(\Delta \mathbf{r})^2>_{||}$ sim', markersize=10)
    # ax.plot(x[:-2], yparam_sim[:-2, 1], 'p',color='r', label=r'$<(\Delta \mathbf{r})^2>_{\perp}$ sim', markersize=10)
    # ax.set_xticks(x[:-2])
    # ax.set_xticklabels(['0', '0.6','0.8', '0.9','0.95'])
    # # ax.legend(loc = 'lower left', fontsize=10)
    # ax.set_ylabel(r'$\alpha^{MC}$', fontsize=fontsize)
    # ax.set_xlabel('Eccentricity', fontsize=fontsize)
    plt.show()
logplot_supp()
# logplot()
# logplot_v5()
################## Linear plot ######################
# fig = plt.figure(figsize=[6,12],tight_layout=True)
# ax1 = fig.add_subplot(311)
# ax2 = fig.add_subplot(312)
# ax3 = fig.add_subplot(313)
# marker = ['o','v','^','<','>','s','+','x']
# def plotter(ax, inmsd_avg):
#     ind = 1
#     for item in inmsd_avg:
#         if ind%2 != 0:
#             x = np.arange(0,len(item))/17
#             ax.plot(x, item/6.25**2, ls='-', marker=marker[ind],label=datalist[ind-1], markevery=2)
#         ind +=1
#     ax.set_xlabel('Lag time (s)', fontsize=15)
#     ax.set_ylabel(r'Average MSD $(\mu m^2)$', fontsize=15)
#     ax.set_xlim([0,50])
#     ax.set_ylim([0,2])
#     ax.legend()
#
# plotter(ax1, msd_avg)
# plotter(ax2, msdx_avg)
# plotter(ax3,msdy_avg)
# plt.show()

