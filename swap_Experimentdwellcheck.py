import os
import matplotlib.pyplot as plt
import module
import numpy as np
import json
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

    n_lines =1 # number of trajs shown in each axis
    if jump == 0:
        dataset = ['ecc0', 'ecc06','ecc08', 'ecc095', 'ecc098','ecc0995']
        ecc = [0,0.6,0.8,0.95,0.98,0.995]
        # dataset = ['ecc0995']
        a = 8 # figure size
        # Create a figure and axis
        fig = plt.figure(figsize=[a,a],tight_layout=True)
        threshold = 0.3 # Pole region is 1/3 of the entire cell length
        ind = 0
        pixelratio = 6.25
        lagtime_lower = 50 # dwell time lower bond in frames. (17fps)
        lagtime_upper = 270 # dwell time upper bond in frames.
        for item in dataset:
            data = tot_vec_overlay_clean[item+'_delx'] # take the major axis position
            datay = tot_vec_overlay_clean[item+'_dely']
            state1 = (data > (np.max(data* threshold)))
            state2 = (data < (np.min(data* threshold)))
            import itertools
            def counts(sig):
                sig = list(sig)
                l = [(k, len(list(g))) for k, g in itertools.groupby(sig)]
                ton = []
                toff = []
                on_ind = [] # start indices of ton
                ind_temp = 0 # current index
                for x in l:
                    if x[0] == 1:
                        ton.append(x[1])
                        on_ind.append(ind_temp)
                    else:
                        toff.append(x[1])
                    ind_temp+=x[1]
                return ton, on_ind

            state1_time, start_ind_s1 = counts(state1)
            state2_time, start_ind_s2 = counts(state2)

            state1_timemask = (np.array(state1_time)<=lagtime_upper)*(np.array(state1_time)>=lagtime_lower) # pick up the indices of the short lag time event
            state1_longtimeind = []
            idx = 0
            for item in state1_timemask:
                if item == True:
                    state1_longtimeind.append(np.arange(start_ind_s1[idx],start_ind_s1[idx]+state1_time[idx]+1))
                idx+=1

            state2_timemask = (np.array(state2_time) <= lagtime_upper)*(np.array(state2_time) >= lagtime_lower)
            state2_longtimeind = []
            idx = 0
            for item in state2_timemask:
                if item == True:
                    state2_longtimeind.append(np.arange(start_ind_s2[idx], start_ind_s2[idx] + state2_time[idx]+1))
                idx += 1

            ax = fig.add_subplot(3,3,ind+1)

            indstart = int(np.random.rand()*len(state1_longtimeind)*0.8)
            for item in state1_longtimeind[indstart:indstart+n_lines]:
                ax.plot(data[item]/pixelratio, datay[item]/pixelratio)
            indstart = int(np.random.rand() * len(state2_longtimeind) * 0.8)
            for item in state2_longtimeind[indstart:indstart+n_lines]:
                ax.plot(data[item]/pixelratio, datay[item]/pixelratio)
            ax.plot(np.max(data* threshold)*np.ones(10)/pixelratio, np.linspace(-2,2,10),'r--')
            ax.plot(np.min(data * threshold)* np.ones(10)/pixelratio, np.linspace(-2, 2, 10), 'r--')
            ax.set_xlim([-2,2])
            ax.set_ylim([-2,2])

            temp_ecc = ecc[ind]
            a,b = module.ellipse_para(temp_ecc)
            theta = np.linspace(0,2*np.pi)
            ax.plot(a*np.cos(theta)*1.05, b*np.sin(theta)*1.05,'k--',label=dataset[ind])
            ax.set_xlabel(r'Position ($\mu m$)')
            ax.set_ylabel(r'Position ($\mu m$)')
            ax.legend()
            ind+=1
        plt.show()
    return
def histlongdwellevent():
    dataset = ['ecc0', 'ecc06','ecc08', 'ecc09', 'ecc095', 'ecc098','ecc0995']
    ecc = [0,0.6,0.8,0.9,0.95,0.98,0.995]
    # dataset = ['ecc0995']
    a = 8 # figure size
    # Create a figure and axis
    fig = plt.figure(figsize=[a,a],tight_layout=True)
    threshold = 0.3 # Pole region is 1/3 of the entire cell length
    ind = 0
    pixelratio = 6.25

    #short dwell event
    dt_s1 = 0 # dwell time lower bond in frames. (17fps)
    dt_s2 = 5 # dwell time upper bond in frames.

    #long dwell event
    lagtime_lower = np.array([1,0.8,0.7,1.1,1.5,2.9,3.3])*17 # Average of long time-scale extracted from the fitting
    lagtime_upper = 270 # dwell time upper bond in frames.

    for item in dataset:
        if item == 'ecc0':
            x = tot_vec_overlay_clean[item + '_delx']  # take the major axis position
            y = tot_vec_overlay_clean[item + '_dely']

            os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/ecc0/')
            a = 0.11 * (np.loadtxt('20210205_data_27.txt')[:, 0] - np.loadtxt('20210205_data_27_ct.txt')[1])
            b = 0.11 * (np.loadtxt('20210205_data_27.txt')[:, 1] - np.loadtxt('20210205_data_27_ct.txt')[0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data

            a = 0.11 * (np.loadtxt('20210205_data_28.txt')[:, 0] - np.loadtxt('20210205_data_28_ct.txt')[1])
            b = 0.11 * (np.loadtxt('20210205_data_28.txt')[:, 1] - np.loadtxt('20210205_data_28_ct.txt')[0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data

            a = 0.11 * (np.loadtxt('20210208_data_2ur.txt')[:, 0] - np.loadtxt('20210208_data_2ur_ct.txt')[
                1])
            b = 0.11 * (np.loadtxt('20210208_data_2ur.txt')[:, 1] - np.loadtxt('20210208_data_2ur_ct.txt')[
                0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data

            a = 0.11 * (np.loadtxt('20210208_data_2bl.txt')[:, 0] - np.loadtxt('20210208_data_2bl_ct.txt')[
                1])
            b = 0.11 * (np.loadtxt('20210208_data_2bl.txt')[:, 1] - np.loadtxt('20210208_data_2bl_ct.txt')[
                0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data
            data = x
            datay = y
        else:
            data = tot_vec_overlay_clean[item + '_delx']  # take the major axis position
            datay = tot_vec_overlay_clean[item + '_dely']
        state1 = (data > (np.max(data* threshold)))
        state2 = (data < (np.min(data* threshold)))
        import itertools
        def counts(sig):
            sig = list(sig)
            l = [(k, len(list(g))) for k, g in itertools.groupby(sig)]
            ton = []
            toff = []
            on_ind = [] # start indices of ton
            ind_temp = 0 # current index
            for x in l:
                if x[0] == 1:
                    ton.append(x[1])
                    on_ind.append(ind_temp)
                else:
                    toff.append(x[1])
                ind_temp+=x[1]
            return ton, on_ind

        state1_time, start_ind_s1 = counts(state1)
        state2_time, start_ind_s2 = counts(state2)

        state1_timemask = (np.array(state1_time)<=lagtime_upper)*(np.array(state1_time)>=lagtime_lower[ind]) # pick up the indices of the short lag time event
        state1_longtimeind = []
        idx = 0
        for item in state1_timemask:
            if item == True:
                state1_longtimeind.append(np.arange(start_ind_s1[idx]+1,start_ind_s1[idx]+state1_time[idx]+1))
            idx+=1

        state2_timemask = (np.array(state2_time) <= lagtime_upper)*(np.array(state2_time) >= lagtime_lower[ind])
        state2_longtimeind = []
        idx = 0
        for item in state2_timemask:
            if item == True:
                state2_longtimeind.append(np.arange(start_ind_s2[idx]+1, start_ind_s2[idx] + state2_time[idx]+1))
            idx += 1

        ax = fig.add_subplot(3,3,ind+1)

        # merge all index into 1D array
        total_longind = []
        for item in state1_longtimeind:
            for number in item:
                total_longind.append(number)

        for item in state2_longtimeind:
            for number in item:
                total_longind.append(number)
        total_longind = np.array(total_longind)-1
        ax.hist2d(data[total_longind]/pixelratio, datay[total_longind]/pixelratio, bins=30, range=[[-2,2], [-2,2]])
        ax.plot(np.max(data* threshold)*np.ones(10)/pixelratio, np.linspace(-2,2,10),'r--')
        ax.plot(np.min(data * threshold)* np.ones(10)/pixelratio, np.linspace(-2, 2, 10), 'r--')
        ax.set_xlim([-2,2])
        ax.set_ylim([-2,2])

        temp_ecc = ecc[ind]
        a,b = module.ellipse_para(temp_ecc)
        theta = np.linspace(0,2*np.pi)
        ax.plot(a*np.cos(theta)*1.05, b*np.sin(theta)*1.05,'k--',label=dataset[ind])
        ax.set_xlabel(r'Position ($\mu m$)',fontsize=15)
        ax.set_ylabel(r'Position ($\mu m$)',fontsize=15)
        ax.legend()
        ind+=1
    plt.show()
def histshortdwellevent():
    dataset = ['ecc0', 'ecc06','ecc08', 'ecc09', 'ecc095', 'ecc098','ecc0995']
    ecc = [0,0.6,0.8,0.9,0.95,0.98,0.995]
    # dataset = ['ecc0995']
    a = 8 # figure size
    # Create a figure and axis
    fig = plt.figure(figsize=[a,a],tight_layout=True)
    threshold = 0.3 # Pole region is 1/3 of the entire cell length
    ind = 0
    pixelratio = 6.25

    #short dwell event
    dt_s1 = 0 # dwell time lower bond in frames. (17fps)
    dt_s2 = 5 # dwell time upper bond in frames.

    #long dwell event
    lagtime_lower = np.array([0.5,0.8,0.7,1.1,1.5,2.9,3.3])*17 # Average of long time-scale extracted from the fitting
    lagtime_upper = 270 # dwell time upper bond in frames.

    for item in dataset:
        if item == 'ecc0':
            x = tot_vec_overlay_clean[item + '_delx']  # take the major axis position
            y = tot_vec_overlay_clean[item + '_dely']

            os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/ecc0/')
            a = 0.11 * (np.loadtxt('20210205_data_27.txt')[:, 0] - np.loadtxt('20210205_data_27_ct.txt')[1])
            b = 0.11 * (np.loadtxt('20210205_data_27.txt')[:, 1] - np.loadtxt('20210205_data_27_ct.txt')[0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data

            a = 0.11 * (np.loadtxt('20210205_data_28.txt')[:, 0] - np.loadtxt('20210205_data_28_ct.txt')[1])
            b = 0.11 * (np.loadtxt('20210205_data_28.txt')[:, 1] - np.loadtxt('20210205_data_28_ct.txt')[0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data

            a = 0.11 * (np.loadtxt('20210208_data_2ur.txt')[:, 0] - np.loadtxt('20210208_data_2ur_ct.txt')[
                1])
            b = 0.11 * (np.loadtxt('20210208_data_2ur.txt')[:, 1] - np.loadtxt('20210208_data_2ur_ct.txt')[
                0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data

            a = 0.11 * (np.loadtxt('20210208_data_2bl.txt')[:, 0] - np.loadtxt('20210208_data_2bl_ct.txt')[
                1])
            b = 0.11 * (np.loadtxt('20210208_data_2bl.txt')[:, 1] - np.loadtxt('20210208_data_2bl_ct.txt')[
                0])
            mask = (np.abs(a) <= 1) * (np.abs(b) <= 1)
            x = np.append(x, 6.25 * a[mask])  # newly added data
            y = np.append(y, 6.25 * b[mask])  # newly added data
            data = x
            datay = y
        else:
            data = tot_vec_overlay_clean[item+'_delx'] # take the major axis position
            datay = tot_vec_overlay_clean[item+'_dely']
        state1 = (data > (np.max(data* threshold)))
        state2 = (data < (np.min(data* threshold)))
        import itertools
        def counts(sig):
            sig = list(sig)
            l = [(k, len(list(g))) for k, g in itertools.groupby(sig)]
            ton = []
            toff = []
            on_ind = [] # start indices of ton
            ind_temp = 0 # current index
            for x in l:
                if x[0] == 1:
                    ton.append(x[1])
                    on_ind.append(ind_temp)
                else:
                    toff.append(x[1])
                ind_temp+=x[1]
            return ton, on_ind

        state1_time, start_ind_s1 = counts(state1)
        state2_time, start_ind_s2 = counts(state2)

        state1_timemask = (np.array(state1_time)<=dt_s2)*(np.array(state1_time)>=dt_s1) # pick up the indices of the short lag time event
        state1_longtimeind = []
        idx = 0
        for item in state1_timemask:
            if item == True:
                state1_longtimeind.append(np.arange(start_ind_s1[idx]+1,start_ind_s1[idx]+state1_time[idx]+1))
            idx+=1

        state2_timemask = (np.array(state2_time) <= dt_s2)*(np.array(state2_time) >= dt_s1)
        state2_longtimeind = []
        idx = 0
        for item in state2_timemask:
            if item == True:
                state2_longtimeind.append(np.arange(start_ind_s2[idx]+1, start_ind_s2[idx] + state2_time[idx]+1))
            idx += 1

        ax = fig.add_subplot(3,3,ind+1)

        # merge all index into 1D array
        total_longind = []
        for item in state1_longtimeind:
            for number in item:
                total_longind.append(number)

        for item in state2_longtimeind:
            for number in item:
                total_longind.append(number)
        total_longind = np.array(total_longind)-1


        ax.hist2d(data[total_longind]/pixelratio, datay[total_longind]/pixelratio, bins=30, range=[[-2,2], [-2,2]])
        ax.plot(np.max(data* threshold)*np.ones(10)/pixelratio, np.linspace(-2,2,10),'r--')
        ax.plot(np.min(data * threshold)* np.ones(10)/pixelratio, np.linspace(-2, 2, 10), 'r--')
        ax.set_xlim([-2,2])
        ax.set_ylim([-2,2])

        temp_ecc = ecc[ind]
        a,b = module.ellipse_para(temp_ecc)
        theta = np.linspace(0,2*np.pi)
        ax.plot(a*np.cos(theta)*1.05, b*np.sin(theta)*1.05,'k--',label=dataset[ind])
        ax.set_xlabel(r'Position ($\mu m$)',fontsize=15)
        ax.set_ylabel(r'Position ($\mu m$)', fontsize=15)
        ax.legend()
        ind+=1
    plt.show()

histlongdwellevent()