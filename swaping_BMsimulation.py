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

main_path = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/"

def swappingplt():
    jump = 0 # jump to zoom plot
    tau_enable = 1
    if jump == 0:
        dataset = ['ecc0', 'ecc06','ecc08', 'ecc09', 'ecc095']
        # dataset = ['ecc0995']
        a = 8 # figure size
        # Create a figure and axis
        fig = plt.figure(figsize=[a/1.6,a],tight_layout=True)

        ind = 0
        bins = [13,12,12,12,12,10,12]
        # bins = [13, 12, 12, 12]
        # bins = 20
        binscum = [300, 300, 300, 300, 300, 300, 300]
        # maxrange = [50, 70, 100, 130, 200, 250, 300]
        maxrange = [183, 256, 367, 476, 733, 917, 1100]
        # maxrange = [183, 256, 733, 1100]
        # maxrange = [1100]
        # maxrange_cum = [300, 300, 300, 300, 300, 300, 300]
        maxrange_cum = 1100

        legend = ['e=0', 'e=0.6','e=0.8', 'e=0.9', 'e=0.95']
        width = np.array(maxrange)/np.array(bins)*0.9 # bar width
        tot_fit = []
        tot_fitcum = []
        threshold = 0.3 # Pole region is 1/3 of the entire cell length
        width_g = 0.15
        tot_cov = []
        tot_covcum = []
        ax = fig.add_subplot(211)
        cmap = cm.get_cmap("inferno").colors
        colorint = 50
        colorshift = 0.8
        cind = 0
        for item in dataset:
            os.chdir(main_path+item)
            data = np.loadtxt('pos.txt')[:, 0]
            data = np.concatenate((data, np.loadtxt('pos2.txt')[:,0]))# take the major axis position
            state1 = (data > (np.max(data* threshold)))
            state2 = (data < (np.min(data* threshold)))
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

            # hist, bin_edge = np.histogram(time, bins=bins[ind], range = [0, maxrange[ind]]) #plotting
            hist, bin_edge = np.histogram(time, bins=bins[ind], range = [0, maxrange[ind]])  # plotting

            # histcum, bin_edgecum = np.histogram(time, bins=binscum[ind], range=[0, maxrange_cum[ind]]) #cumulative fitting
            histcum, bin_edgecum = np.histogram(time, bins=binscum[ind], range=[0, maxrange_cum])  # cumulative fitting
            cumutive = np.cumsum(histcum)
            # if ind == 0 or ind == 1 or ind == 4 or ind == 6:
            print('ind is:'+str(ind))
            ax.bar(bin_edge[:-1], hist/hist[0],color=cmap[int((cind+colorshift)*colorint)],align = 'edge', width=width[ind], edgecolor = 'k', zorder=-ind,alpha=0.8)
            # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/bin', bin_edge)
            # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/hist',  hist/hist[0])
            # print('a')
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
                                  p0=np.array([1, 1, 10, 200,1]),
                                  bounds=(0,np.inf))
            print('poptcum:'+str(poptcum))
            tot_fitcum.append(poptcum)
            tot_covcum.append(pcovcum)

            # plot exponential fitting
            xt = np.linspace(0, fitx_x[-1],1000)
            yt = exp_func(xt, popt[0], popt[1], popt[2], popt[3])
            ytcum = exp_func(xt, poptcum[0], poptcum[1], poptcum[2], poptcum[3])
            # if ind == 0 or ind == 1 or ind == 4 or ind == 6:
            ax.plot(xt, ytcum / ytcum[0], '--', zorder=ind, linewidth=2, label=legend[cind],
                    color=cmap[int((cind + colorshift-0.2) * colorint)])
            cind += 1
            if item =='ecc095':
                ax.set_yticks([10, 100, 1000])
            # xtick = np.linspace(0, maxrange[-1]-maxrange[-1]*0.1,3)
            # xticklabel = np.around(xtick, decimals=1)
            # ax.set_xticks(xtick)
            # ax.set_xticklabels(xticklabel)
            ax.set_yscale('log')
            ax.set_ylim([1e-4, 1.3])
            ax.set_xlabel('Dwell time', fontsize=15)
            ax.set_ylabel('Normalized frequency', fontsize=15)
            ax.legend()
            ax.set_xlim([-2,800])
            ind += 1

        ############
        # tau plot #
        ############
        if tau_enable == 1:
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
            # x = [0]
            # y = [0, 2, 4]
            ax.errorbar(x[:-2], ts, yerr=es, linestyle='None', marker='.', ms=13, capsize=8, color='r', label='Short dwell time')
            ax.errorbar(x[:-2], tl, yerr=el, linestyle='None', marker='.', ms=13, capsize=8, color='k', label='Long dwell time')
            np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/ts',ts)
            np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/es', es)
            np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/tl', tl)
            np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/el', el)
            ax.set_xlabel('Cavity eccentricity', fontsize=15)
            ax.set_ylabel('Dwell time', fontsize=15)
            ax.set_xticks([0,2,4], minor=False)
            # ax.set_xticks([1, 2, 4,5], minor=True)

            # ax.set_yticks(y)
            ax.set_xticklabels(['0', '0.8', '0.95'])
            ax.tick_params(axis='both', which='both', labelsize=ticklabelsize)
            ax.legend(fontsize=15)
        plt.show()


    # else:
    #     ##########Zoom in plot #####
    #     fig = plt.figure()
    #     maxrange = 30
    #     bins = 13
    #     threshold = 0.3
    #     ax = fig.add_subplot(111)
    #     data = tot_vector['ecc0995_delx']
    #     state1 = (data>(np.max(data)*threshold))
    #     state2 = (data <(np.min(data) * threshold))
    #
    #     import itertools
    #     def counts(sig):
    #         sig = list(sig)
    #         l = [(k, len(list(g))) for k, g in itertools.groupby(sig)]
    #         ton = []
    #         toff = []
    #         for x in l:
    #             if x[0] == 1:
    #                 ton.append(x[1])
    #             else:
    #                 toff.append(x[1])
    #         return ton
    #
    #     state1_time = counts(state1)
    #     state2_time = counts(state2)
    #     time = state1_time+state2_time
    #     time = np.array(time)
    #     # time = time[time>1]
    #     hist, bin_edge = np.histogram(time, bins=bins, range = [0, maxrange])
    #     ax.bar(bin_edge[:-1], hist, align = 'center', width=maxrange/bins*0.9, color = 'r', edgecolor = 'k', label="Ecc=0.995")
    #     ax.set_yscale('log')
    #
    #     xtick = np.linspace(0, maxrange - maxrange * 0.1, 3)
    #     xticklabel = np.around(xtick / 17, decimals=1)
    #     ax.set_xticks(xtick)
    #     ax.set_xticklabels(xticklabel)
    #     ax.set_xlabel('Dwell time (s)', fontsize=16)
    #     ax.set_ylabel('Counts', fontsize=16)
    #     ax.tick_params(axis='both', labelsize=15)
    #
    #     plt.show()
    #     # Exp fitting
    #     # Linear-exponential fit method
    #     # fitx_x = bin_edge[start:-1]
    #     # hist1 = hist[start:]
    #     # mask = (hist1!=0)
    #     # fit_para, cov = np.polyfit(fitx_x[mask], np.log(hist1[mask]), deg=int(1), full=False, cov=True)
    #     # tot_fit.append(fit_para)
    #     # tot_cov.append(cov)
    return
def plotter():
    x = np.loadtxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/ecc0995/pos_polejump.txt')
    plt.hist2d(x[:,0], x[:,1], bins=30)
    # plt.plot(x[:1000,0],x[:1000,1])
    plt.show()
swappingplt()