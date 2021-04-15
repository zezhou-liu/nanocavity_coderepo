import tifffile as tif
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
def cavcenterfromavg():
    # Directory setup
    data_savepth = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/'
    savefolder = '/ecc0'
    vidpath = '/media/zezhou/Elements/DATA_McGill/2021_Feb/20210208/data_16/'
    avgpath = 'AVG_Substack (1-7000)_background-subtracted.tif'
    centerfile = '20210208_data_16_ct.txt'
    avg = tif.imread(vidpath+avgpath)
    os.chdir(data_savepth+savefolder)
    threshold = 0.02 # pixel brighter than threshold*maxintensity will be counted.
    mask = avg>threshold*np.max(avg)

    xx, yy = np.meshgrid(range(mask.shape[1]), range(mask.shape[0]))

    cx = np.sum(np.multiply(xx, mask))/np.sum(mask)
    cy = np.sum(np.multiply(yy, mask))/np.sum(mask)
    center = np.array([cx, cy])
    np.savetxt(centerfile, center)


    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.pcolormesh(xx, yy,avg)
    ax.plot(cx,cy,'x')
    ax = fig.add_subplot(212)
    ax.pcolormesh(xx, yy,mask)
    ax.plot(cx,cy,'x')
    plt.show()
    print('done')

data_savepth = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/'
savefolder = '/ecc0'
os.chdir(data_savepth+savefolder)

# '20210205_data_23',
filelist = ['20210205_data_27', '20210205_data_28',
            '20210208_data_2ur', '20210208_data_2bl']
# 0.042,
# tof = [0.0207651,0.0210107,0.0210929,0.0210929]#time of frame in sec
tof = 0.021 # in sec
maxt = 100 # max time in sec
pixelsize = 0.11 # size of a single pixel. 0.11 um in this case.
time = 1 # fitting time in sec
def msdcal(filelist):
    MSD = []
    MSDx = []
    MSDy = []
    for item in filelist:
        filename = item+'.txt'
        centername = item+'_ct.txt'
        data = np.loadtxt(filename)

        for i in range(data.shape[0]):
            if data[i,0]==0:
                data[i,:]=(data[i-1,:]+data[i+1,:])/2

        center = np.loadtxt(centername)

        data[:,0] = data[:,0]-center[1]
        data[:,1] = data[:,1]-center[0]
        xtemp = data[:,0]
        ytemp = data[:,1]
        l = xtemp.size
        deno = np.arange(1,l+1)
        deno = np.concatenate((deno, np.flip(deno[0:-1])))
        MSD.append(2*(np.mean(xtemp**2+ytemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno)-2*np.divide(np.correlate(ytemp,ytemp,mode='full'), deno))
        MSDx.append(2*(np.mean(xtemp**2))-2*np.divide(np.correlate(xtemp,xtemp,mode='full'), deno))
        MSDy.append(2 * (np.mean(ytemp ** 2)) - 2 * np.divide(np.correlate(ytemp, ytemp, mode='full'), deno))
    return MSD, MSDx, MSDy
def msd_average(inmsd):
    nframes = []
    nclip=len(inmsd)
    for index in range(nclip):
        nframes.append(inmsd[index].size)
    nframes_in = nframes
    arr = np.zeros((len(nframes_in),
                    np.max(nframes_in)))  # note msd is symmetric and we only process the positive(right side) value.
    for idx in range(len(inmsd)):
        inmsd[idx][[int((len(inmsd[idx]) - 1) / 2)]] = 1e-15
        arr[idx, int((np.max(nframes_in) - 1) / 2):int((np.max(nframes_in) - 1) / 2 + (nframes_in[idx] + 1) / 2)] = \
        inmsd[idx][int((nframes_in[idx] - 1) / 2):]  # save the right side value.

    # sort array from min length to max length
    indarray = np.argsort(nframes_in)
    nframes_in = np.array(nframes_in)[indarray]
    arr = arr[indarray, :]
    # weight list
    msd_temp = np.zeros(int((np.max(nframes_in) + 1) / 2))  # mean value of msd
    std_temp = np.zeros(int((np.max(nframes_in) + 1) / 2))  # std of msd

    stdarr = arr
    stdarr[stdarr == 0] = np.nan
    for idx in range(len(inmsd)):
        if idx == 0:
            msd_temp[0:int((nframes_in[idx] + 1) / 2)] = np.average(
                arr[:, int((np.max(nframes_in) - 1) / 2):int((np.max(nframes_in) - 1) / 2 + (nframes_in[idx] + 1) / 2)],
                axis=0, weights=nframes_in)
            std_temp[0:int((nframes_in[idx] + 1) / 2)] = np.nanstd(stdarr[:, int((np.max(nframes_in) - 1) / 2):int(
                (np.max(nframes_in) - 1) / 2 + (nframes_in[idx] + 1) / 2)], axis=0)
        else:
            s = int((nframes_in[idx - 1] + 1) / 2)
            e = int((nframes_in[idx] + 1) / 2)
            msd_temp[int((nframes_in[idx - 1] + 1) / 2):int((nframes_in[idx] + 1) / 2)] = np.average(arr[idx:,
                                                                                                     int((np.max(
                                                                                                         nframes_in) - 1) / 2 + (
                                                                                                                     nframes_in[
                                                                                                                         idx - 1] + 1) / 2):
                                                                                                     int((np.max(
                                                                                                         nframes_in) - 1) / 2 + (
                                                                                                                     nframes_in[
                                                                                                                         idx] + 1) / 2)],
                                                                                                     axis=0,
                                                                                                     weights=nframes_in[
                                                                                                             idx:])
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

msd, msdx, msdy = msdcal(filelist)
mean = msd_average(msd)
xmean = msd_average(msdx)
ymean = msd_average(msdy)

np.savetxt('msd_mean', mean)
np.savetxt('msdx_mean', xmean)
np.savetxt('msdy_mean', ymean)
def exponentfit_individualvid(save=False):
    # Diffusion exponent fitting based on each video. This way is supposed to extract the errorbar.
      # fitting time (sec)
    def fitfunc(x, a, b):  # fitting function for MSD. Exponent
        return a * x ** b
    sol = []
    cov = []
    solx = []
    covx = []
    soly = []
    covy = []
    for clip_ind in range(len(msd)):
        ftime = int(time /tof)
        msd_tmp = msd[clip_ind]
        start = int((len(msd_tmp)-1)/2) # take the right side of MSD since it's symmetric
        end = start + ftime
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
        popt, pcov = curve_fit(fitfunc, t_tmp, msd_tmp*pixelsize**2, p0=[1, 1])
        poptx, pcovx = curve_fit(fitfunc, t_tmp, msdx_tmp*pixelsize**2, p0=[1, 1])
        popty, pcovy = curve_fit(fitfunc, t_tmp, msdy_tmp*pixelsize**2, p0=[1, 1])

        sol.append(popt)
        cov.append(pcov)
        solx.append(poptx)
        covx.append(pcovx)
        soly.append(popty)
        covy.append(pcovy)
    sol.append([0,0.683])
    solx.append([0, 0.641])
    soly.append([0, 0.728])
    if save==True:
        os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/ecc0/')
        f = open('msdfit_clip.txt', 'a')
        np.savetxt(f, sol)
        f.close()
        f = open('msdfit_cov_clip.txt', 'a')
        for item11 in cov:
            np.savetxt(f, item11)
        f.close()

        f = open('msdfitx_clip.txt', 'a')
        np.savetxt(f, solx)
        f.close()
        f = open('msdfitx_cov_clip.txt', 'a')
        for item11 in covx:
            np.savetxt(f, item11)
        f.close()

        f = open('msdfity_clip.txt', 'a')
        np.savetxt(f, soly)
        f.close()
        f = open('msdfity_cov_clip.txt', 'a')
        for item11 in covy:
            np.savetxt(f, item11)
        f.close()
    return np.array(sol), np.array(solx), np.array(soly)

sol, solx, soly = exponentfit_individualvid(save=False)

fig = plt.figure(figsize=[5,5])
ax = fig.add_subplot(221)
axx = fig.add_subplot(222)
axy = fig.add_subplot(223)
test = fig.add_subplot(224)

test.plot(xmean[0][:int(time/tof)],'+')


def fitfuncc(x, a, b):  # fitting function for MSD. Exponent
    return a * x ** b

for i in range(len(msd)):
    t = np.arange(msd[i].size)*tof
    start = int((msd[i].size-1)/2)
    end = start + int(maxt/tof)
    a = int(maxt/tof)
    x = t[:int(maxt/tof)]
    print('a')
    ax.plot(t[:int(maxt/tof)], msd[i][start:end]*pixelsize**2,'+')
    axx.plot(t[:int(maxt/tof)], msdx[i][start:end]*pixelsize**2,'+')
    axy.plot(t[:int(maxt/tof)], msdy[i][start:end]*pixelsize**2,'+')

ttemp = np.linspace(0,time,100)
print('Full-Mean:'+str(np.mean(sol[:,1])))
print('SEM:'+str(np.std(sol[:,1])/np.sqrt(len(filelist))))
print('X-Mean:'+str(np.mean(solx[:,1])))
print('SEM:'+str(np.std(solx[:,1])/np.sqrt(len(filelist))))
print('Y-Mean:'+str(np.mean(soly[:,1])))
print('Y SEM:'+str(np.std(solx[:,1])/np.sqrt(len(filelist))))

ax.plot(ttemp, fitfuncc(ttemp, np.mean(sol[:,0]), np.mean(sol[:,1])))
axx.plot(ttemp, fitfuncc(ttemp, np.mean(solx[:,0]), np.mean(solx[:,1])))
axy.plot(ttemp, fitfuncc(ttemp, np.mean(soly[:,0]), np.mean(soly[:,1])))

# ax.set_xscale('log')
# ax.set_yscale('log')
# axx.set_xscale('log')
# axx.set_yscale('log')
# axy.set_xscale('log')
# axy.set_yscale('log')
# ax.set_xlim([6e-2, 100])
# axx.set_xlim([6e-2, 100])
# axy.set_xlim([6e-2, 100])
# ax.set_ylim([6e-2, 100])
# axx.set_ylim([6e-2, 100])
# axy.set_ylim([6e-2, 100])
plt.show()