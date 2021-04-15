import numpy as np
import matplotlib.pyplot as plt
import module
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import Normalize
import json
import os
import matplotlib.colors as mcolors
from matplotlib.transforms import Affine2D
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

main_path = "/media/zezhou/Seagate Expansion Drive/McGillResearch/2019Manuscript_Analysis/Analysis/datafterlinearshift/tplasmid/"
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
    handle1, tot_file_shift = module.allshift(handle1) # Shift data to zero according to YOYO-3 channel
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

def swappingplt():
    jump = 0 # jump to zoom plot
    tot_vector = handle1.tot_vec_overlay_clean
    if jump == 0:
        dataset = ['ecc0', 'ecc06','ecc08','ecc09','ecc095', 'ecc098', 'ecc0995']
        a = 13 # figure size
        # Create a figure and axis
        fig = plt.figure(figsize=[a, a/1.3333333])

        ind = 0
        # bins = [13,12,12,12,12,12,12]
        bins = [200, 200, 200, 200, 200, 200, 200]
        # maxrange = [50, 70, 100, 130, 200, 250, 300]
        maxrange = [300, 300, 300, 300, 300, 300, 300]
        pltrange = [1,1,1,1,1,1,1]
        start = 3 # starting bins
        maxtime = np.array(maxrange)/17
        legend = ['Ecc=0', 'Ecc=0.6', 'Ecc=0.8', 'Ecc=0.9', 'Ecc=0.95', 'Ecc=0.98', 'Ecc=0.995']
        width = np.array(maxrange)/np.array(bins)*0.9 # bar width
        tot_fit = []
        # threshold = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        threshold = 0.3*np.ones(7)
        width_g = 0.15
        tot_cov = []
        tot_fit2 = []
        tot_cov2 = []
        # ax = fig.add_subplot(121)
        for item in dataset:
            r = (ind)//3
            c = (ind)%3
            ax = fig.add_axes([0.1+c*(width_g+0.07), 0.7-r*(width_g+0.07), width_g, width_g])
            data = tot_vector[item+'_delx']
            state1 = (data>(np.max(data)*threshold[ind]))
            state2 = (data <(np.min(data) * threshold[ind]))

            import itertools
            def counts(sig):
                sig = list(sig)
                l = [(k, len(list(g))) for k, g in itertools.groupby(sig)]
                ton = []
                toff = []
                for x in l:
                    if x[0] == 1:
                        ton.append(x[1])
                    else:
                        toff.append(x[1])
                return ton

            state1_time = counts(state1)
            state2_time = counts(state2)
            time = np.array(state1_time+state2_time)
            hist, bin_edge = np.histogram(time, bins=bins[ind], range = [0, maxrange[ind]])
            cumutive = np.cumsum(hist)
            ax.plot(bin_edge[:-1], cumutive/cumutive[-1], 'r+', markersize=10, label=legend[ind])
            np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/bin', bin_edge)
            np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/cum', cumutive/cumutive[-1])
            from scipy.optimize import curve_fit
            def cum_func(x,a,b,c,d,e):
                return -a*b*np.exp(-x/b)-c*d*np.exp(-x/d)+e

            # cum fit
            fitx_x = bin_edge[:-1]
            fitx_y = cumutive/cumutive[-1]
            popt,pcov = curve_fit(cum_func, fitx_x, fitx_y,
                                  p0=np.array([1, 10, 1, 10,100]),
                                  bounds=(0,np.inf))

            tot_fit.append(popt)
            print(item+' cumulative solution:'+str(popt))
            tot_cov.append(pcov)
            print(item + ' cumulative covariance:' + str(pcov))
            # plot exponential fitting
            xt = np.linspace(0, fitx_x[-1],1000)
            ax.plot(xt[:int(pltrange[ind]*1000)], cum_func(xt[:int(pltrange[ind]*1000)], popt[0], popt[1], popt[2], popt[3],popt[4]), 'k--')
            if item =='ecc095':
                ax.set_yticks([10, 100, 1000])
            # ax.legend()
            xtick = np.linspace(0, maxrange[ind]-maxrange[ind]*0.1,3)
            xticklabel = np.around(xtick/17, decimals=1)
            ax.set_xticks(xtick)
            ax.set_xticklabels(xticklabel)
            ax.set_yscale('log')
            ax.set_xlabel('Dwell time (s)', fontsize=13)
            ax.set_ylabel('Counts', fontsize=13)
            ax.legend(fontsize=12, loc='lower right')
            ind += 1

        ############
        # tau plot #
        ############
        # fig = plt.figure(figsize=[6.4, 6.4 / 1.333])
        # ax = fig.add_axes([0.2, 0.2, 0.5, 0.5])
        r = 7 // 3
        c = 7 % 3
        ax = fig.add_axes([0.1+c*(width_g+0.07), 0.7-r*(width_g+0.07), width_g, width_g])
        # ax = fig.add_subplot(122)
        labelsize = 13
        ticklabelsize = 11
        time = []
        error = []

        time1 = []
        error1 = []

        time2 = []
        error2 = []
        for i in range(len(tot_fit)):
            # time.append(-1 / tot_fit[i][0] / 17)
            # error.append(1 / (tot_fit[i][0] ** 2) * np.sqrt(tot_cov[i][0][0]) / 17)
            time1.append(tot_fit[i][1] / 17)
            error1.append(np.sqrt(tot_cov[i][1][1]) / 17)

            time2.append(tot_fit[i][3] / 17)
            error2.append(np.sqrt(tot_cov[i][3][3]) / 17)
        ts = np.minimum(time1, time2)
        tl = np.maximum(time1, time2)
        es = []
        el = []
        inda = len(time1)
        for i in range(inda):
            if time1[i]>time2[i]:
                el.append(error1[i])
                es.append(error2[i])
            else:
                el.append(error2[i])
                es.append(error1[i])
        x = [0, 1, 2, 3,4,5,6]
        # y = [0, 2, 4]
        ax.errorbar(x, ts, yerr=es, linestyle='None', marker='.', ms=13, capsize=8, color='r', label='Short dwell time')
        ax.errorbar(x, tl, yerr=el, linestyle='None', marker='.', ms=13, capsize=8, color='k', label='Long dwell time')
        # ax.plot(x,ts,'ro')
        # ax.plot(x, tl, 'ko')
        ax.set_xlabel('Cavity eccentricity', fontsize=labelsize)
        ax.set_ylabel('Dwell time (s)', fontsize=labelsize)
        ax.set_xticks([0,3,6], minor=False)
        ax.set_xticks([1, 2, 4,5], minor=True)
        # ax.set_yticks(y)
        ax.set_xticklabels(['0', '0.9', '0.995'])
        ax.tick_params(axis='both', which='both', labelsize=ticklabelsize)
        ax.legend(fontsize=10, loc='upper left', borderpad=0.1)
        r = 8 // 3
        c = 8 % 3
        ax = fig.add_axes([0.1 + c * (width_g + 0.07), 0.7 - r * (width_g + 0.07), width_g, width_g])
        # ax = fig.add_subplot(122)
        labelsize = 13
        ticklabelsize = 11
        time = []
        error = []

        time1 = []
        error1 = []

        time2 = []
        error2 = []
        for i in range(len(tot_fit)):
            # time.append(-1 / tot_fit[i][0] / 17)
            # error.append(1 / (tot_fit[i][0] ** 2) * np.sqrt(tot_cov[i][0][0]) / 17)
            time1.append(tot_fit[i][1] / 17)
            error1.append(np.sqrt(tot_cov[i][1][1]) / 17)

            time2.append(tot_fit[i][3] / 17)
            error2.append(np.sqrt(tot_cov[i][3][3]) / 17)
        ts = np.minimum(time1, time2)
        tl = np.maximum(time1, time2)
        es = []
        el = []
        inda = len(time1)
        for i in range(inda):
            if time1[i] > time2[i]:
                el.append(error1[i])
                es.append(error2[i])
            else:
                el.append(error2[i])
                es.append(error1[i])
        print('time long:'+str(tl))
        print('time short:' + str(ts))
        x = [0, 1, 2, 3, 4, 5, 6]
        # y = [0, 2, 4]
        ax.errorbar(x, ts, yerr=es, linestyle='None', marker='.', ms=13, capsize=8, color='r')
        ax.errorbar(x, tl, yerr=el, linestyle='None', marker='.', ms=13, capsize=8, color='k')
        ax.set_xlabel('Cavity eccentricity', fontsize=labelsize)
        ax.set_ylabel('Dwell time (s)', fontsize=labelsize)
        ax.set_xticks([0, 3, 6], minor=False)
        ax.set_xticks([1, 2, 4, 5], minor=True)
        # ax.set_yticks(y)
        ax.set_xticklabels(['0', '0.9', '0.995'])
        ax.tick_params(axis='both', which='both', labelsize=ticklabelsize)
        ax.set_ylim(0,0.5)

        plt.show()


    else:
        ##########Zoom in plot #####
        fig = plt.figure()
        maxrange = 30
        bins = 13
        threshold = 0.3
        ax = fig.add_subplot(111)
        data = tot_vector['ecc0995_delx']
        state1 = (data>(np.max(data)*threshold))
        state2 = (data <(np.min(data) * threshold))

        import itertools
        def counts(sig):
            sig = list(sig)
            l = [(k, len(list(g))) for k, g in itertools.groupby(sig)]
            ton = []
            toff = []
            for x in l:
                if x[0] == 1:
                    ton.append(x[1])
                else:
                    toff.append(x[1])
            return ton

        state1_time = counts(state1)
        state2_time = counts(state2)
        time = state1_time+state2_time
        time = np.array(time)
        # time = time[time>1]
        hist, bin_edge = np.histogram(time, bins=bins, range = [0, maxrange])
        ax.bar(bin_edge[:-1], hist, align = 'center', width=maxrange/bins*0.9, color = 'r', edgecolor = 'k', label="Ecc=0.995")
        ax.set_yscale('log')

        xtick = np.linspace(0, maxrange - maxrange * 0.1, 3)
        xticklabel = np.around(xtick / 17, decimals=1)
        ax.set_xticks(xtick)
        ax.set_xticklabels(xticklabel)
        ax.set_xlabel('Dwell time (s)', fontsize=16)
        ax.set_ylabel('Counts', fontsize=16)
        ax.tick_params(axis='both', labelsize=15)

        plt.show()
        # Exp fitting
        # Linear-exponential fit method
        # fitx_x = bin_edge[start:-1]
        # hist1 = hist[start:]
        # mask = (hist1!=0)
        # fit_para, cov = np.polyfit(fitx_x[mask], np.log(hist1[mask]), deg=int(1), full=False, cov=True)
        # tot_fit.append(fit_para)
        # tot_cov.append(cov)
    return
swappingplt()