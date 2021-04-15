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
        a = 8 # figure size
        # Create a figure and axis
        fig = plt.figure(figsize=[a/1.6,a],tight_layout=True)

        ind = 0
        bins = [12,12,12,12,12,10,12]
        binscum = [200, 200, 200, 200, 200, 200, 200]
        maxrange = [50, 70, 100, 130, 200, 250, 300]
        maxrange_cum = [300, 300, 300, 300, 300, 300, 300]
        maxtime = np.array(maxrange)/17
        legend = ['e=0', 'e=0.6', 'e=0.8', 'e=0.9', 'e=0.95', 'e=0.98', 'e=0.995']
        width = np.array(maxrange)/np.array(bins)*0.9 # bar width
        tot_fit = []
        tot_fitcum = []
        threshold = 0.3*np.ones(7)
        width_g = 0.15
        tot_cov = []
        tot_covcum = []
        ax = fig.add_subplot(211)
        cmap = cm.get_cmap("inferno").colors
        colorint = 62
        colorshift = 0.8
        cind=0
        for item in dataset:
            if item == 'ecc0':
                data = tot_vector[item+'_delx']
                # frame: 3388
                # newly added data
                os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/ecc0/')
                a = 0.11 * (np.loadtxt('20210205_data_27.txt')[:, 0] - np.loadtxt('20210205_data_27_ct.txt')[1])
                b = 0.11 * (np.loadtxt('20210205_data_27.txt')[:, 1] - np.loadtxt('20210205_data_27_ct.txt')[0])
                mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
                data = np.append(data, 6.25 * a[mask])  # newly added data

                a = 0.11 * (np.loadtxt('20210205_data_28.txt')[:, 0] - np.loadtxt('20210205_data_28_ct.txt')[1])
                b = 0.11 * (np.loadtxt('20210205_data_28.txt')[:, 1] - np.loadtxt('20210205_data_28_ct.txt')[0])
                mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
                data = np.append(data, 6.25 * a[mask])  # newly added data

                a = 0.11 * (np.loadtxt('20210208_data_2ur.txt')[:, 0] - np.loadtxt('20210208_data_2ur_ct.txt')[
                    1])
                b = 0.11 * (np.loadtxt('20210208_data_2ur.txt')[:, 1] - np.loadtxt('20210208_data_2ur_ct.txt')[
                    0])
                mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
                data = np.append(data, 6.25 * a[mask])  # newly added data

                a = 0.11 * (np.loadtxt('20210208_data_2bl.txt')[:, 0] - np.loadtxt('20210208_data_2bl_ct.txt')[
                    1])
                b = 0.11 * (np.loadtxt('20210208_data_2bl.txt')[:, 1] - np.loadtxt('20210208_data_2bl_ct.txt')[
                    0])
                mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
                data = np.append(data, 6.25 * a[mask])  # newly added data
            else:
                data = tot_vector[item+'_delx']

            if item == 'ecc0':
                f1 = 3388 # number of frames from the old data
                state1 = (data[:f1] > (np.max(data[:f1]) * threshold[ind]))
                state2 = (data[:f1] < (np.min(data[:f1]) * threshold[ind]))

                state1 = np.append(state1, data[f1:] > (np.max(data[f1:]) * threshold[ind])) # newly added data. frame rate is different.
                state2 = np.append(state2, data[f1:] < (np.min(data[f1:]) * threshold[ind]))
            else:
                state1 = (data > (np.max(data) * threshold[ind]))
                state2 = (data < (np.min(data) * threshold[ind]))
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
                return np.array(ton)
            if item == 'ecc0':
                f1 = 3388 # number of frames from the old data
                state1_time = counts(state1[:f1]) # old data
                state2_time = counts(state2[:f1])

                state1_time = np.append(state1_time, counts(state1[f1:])*0.021*17) # new data
                state2_time = np.append(state1_time, counts(state2[f1:])*0.021*17)
            else:
                state1_time = counts(state1)
                state2_time = counts(state2)
            time = np.array(state1_time.tolist()+state2_time.tolist())
            # time = time[time>1] # delete 1 frame event?
            hist, bin_edge = np.histogram(time, bins=bins[ind], range = [0, maxrange[ind]]) #plotting

            histcum, bin_edgecum = np.histogram(time, bins=binscum[ind], range=[0, maxrange_cum[ind]]) #cumulative fitting

            cumutive = np.cumsum(histcum)
            if ind == 0 or ind == 1 or ind == 4 or ind == 6:
                ax.bar(bin_edge[:-1], hist/hist[0],color=cmap[int((cind+colorshift)*colorint)],align = 'edge', width=width[ind], edgecolor = 'k', zorder=-ind,alpha=0.8)
                # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/bin_edge', bin_edge)
                # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/hist',
                #            hist/hist[0])
                print('1')
            # Curve_fit method
            from scipy.optimize import curve_fit
            def exp_func(x, a, b, c, d):
                return a*np.exp(-x/b)+c*np.exp(-x/d)
            def cum_func(x,a,b,c,d,e):
                return -a*b*np.exp(-x/b)-c*d*np.exp(-x/d)+e
            fitx_x = bin_edge[:-1]
            fitx_y = hist

            popt, pcov = curve_fit(exp_func, fitx_x, fitx_y,
                                   p0=np.array([1,1,1,1]),
                                   bounds=(0,np.inf))
            print(str(ind)+' popt:'+str(popt))
            # print(pcov)
            # print('check 1')
            tot_fit.append(popt)
            tot_cov.append(pcov)

            # cum fit
            fitx_x = bin_edgecum[:-1]
            fitx_y = cumutive/cumutive[-1]
            poptcum,pcovcum = curve_fit(cum_func, fitx_x, fitx_y,
                                  p0=np.array([1, 1, 1, 1,100]),
                                  bounds=(0,np.inf))
            print('poptcum:'+str(poptcum))
            tot_fitcum.append(poptcum)
            tot_covcum.append(pcovcum)

            # plot exponential fitting
            xt = np.linspace(0, fitx_x[-1],1000)
            yt = exp_func(xt, popt[0], popt[1], popt[2], popt[3])
            ytcum = exp_func(xt, poptcum[0], poptcum[1], poptcum[2], poptcum[3])
            if ind == 0 or ind == 1 or ind ==4 or ind == 6:
                ax.plot(xt, ytcum / ytcum[0], '--', zorder=ind, linewidth=2, label=legend[ind],
                        color=cmap[int((cind + colorshift-0.2) * colorint)])
                cind+=1
            if item =='ecc095':
                ax.set_yticks([10, 100, 1000])
            # ax.legend()
            xtick = np.linspace(0, maxrange[ind]-maxrange[ind]*0.1,3)
            xticklabel = np.around(xtick/17, decimals=1)
            ax.set_xticks(xtick)
            ax.set_xticklabels(xticklabel)
            ax.set_yscale('log')
            ax.set_xlabel('Dwell time (s)', fontsize=15)
            ax.set_ylabel('Normalized frequency', fontsize=15)
            ax.legend()
            ind += 1

        ############
        # tau plot #
        ############
        # fig = plt.figure(figsize=[6.4, 6.4 / 1.333])
        # ax = fig.add_axes([0.2, 0.2, 0.5, 0.5])
        # r = 7 // 3
        # c = 7 % 3
        # ax = fig.add_axes([0.1+c*(width_g+0.07), 0.7-r*(width_g+0.07), width_g, width_g])
        ax = fig.add_subplot(212)
        labelsize = 13
        ticklabelsize = 11
        time = []
        error = []

        time1 = []
        error1 = []

        time2 = []
        error2 = []

        # for i in range(len(tot_fit)):
        #     time.append(-1 / tot_fit[i][0] / 17)
        #     error.append(1 / (tot_fit[i][0] ** 2) * np.sqrt(tot_cov[i][0][0]) / 17)
        #     time1.append(tot_fit[i][1] / 17)
        #     error1.append(np.sqrt(tot_cov[i][1][1]) / 17)
        #
        #     time2.append(tot_fit[i][3] / 17)
        #     error2.append(np.sqrt(tot_cov[i][3][3]) / 17)

        for i in range(len(tot_fitcum)):
            time.append(-1 / tot_fitcum[i][0] / 17)
            error.append(1 / (tot_fitcum[i][0] ** 2) * np.sqrt(tot_covcum[i][0][0]) / 17)
            time1.append(tot_fitcum[i][1] / 17)
            error1.append(np.sqrt(tot_covcum[i][1][1]) / 17)

            time2.append(tot_fitcum[i][3] / 17)
            error2.append(np.sqrt(tot_covcum[i][3][3]) / 17)
            # time1.append(tot_fitcum[i][1] / 17)
            # error1.append(np.sqrt(tot_covcum[i][1][1]) / 17)
            #
            # time2.append(tot_fitcum[i][3] / 17)
            # error2.append(np.sqrt(tot_covcum[i][3][3]) / 17)
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
        print('time long:' + str(tl))
        print('time short:' + str(ts))
        x = [0, 1, 2, 3,4,5,6]
        # y = [0, 2, 4]
        ax.errorbar(x, ts, yerr=es, linestyle='None', marker='.', ms=13, capsize=8, color='r', label='Short dwell time')
        ax.errorbar(x, tl, yerr=el, linestyle='None', marker='.', ms=13, capsize=8, color='k', label='Long dwell time')

        np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/ts', ts)
        np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/es', es)
        np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/tl', tl)
        np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/el', el)
        ax.set_xlabel('Cavity eccentricity', fontsize=15)
        ax.set_ylabel('Dwell time (s)', fontsize=15)
        ax.set_xticks([0,3,6], minor=False)
        ax.set_xticks([1, 2, 4,5], minor=True)
        # ax.set_yticks(y)
        ax.set_xticklabels(['0', '0.9', '0.995'])
        ax.tick_params(axis='both', which='both', labelsize=ticklabelsize)
        ax.legend(fontsize=15)
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