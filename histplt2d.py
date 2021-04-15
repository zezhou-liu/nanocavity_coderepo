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

################################## Visualization #########################################

###############
# histogram2d #
###############
# Details of the plot data are in function. The unlisted data are directly from datahandle.
##########################
# T4-Plasmid
def ecc0_tp():
    x = np.array([])
    y = np.array([])
    x = np.append(x, handle1.tot_file_shift['ecc0_1_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc0_1_y1y'])
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/additionaldata/ecc0/')
    a = 0.11*(np.loadtxt('20210205_data_27.txt')[:, 0] - np.loadtxt('20210205_data_27_ct.txt')[1])
    b = 0.11*(np.loadtxt('20210205_data_27.txt')[:, 1] - np.loadtxt('20210205_data_27_ct.txt')[0])
    mask = (np.abs(a)<=1)*(np.abs(b)<=1)
    x = np.append(x, 6.25*a[mask]) # newly added data
    y = np.append(y, 6.25*b[mask])  # newly added data

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
    # x = x - np.mean(x)
    print('tot_pos_overlay_shift ecc0:' + str(np.mean(x)))
    return x, y
def ecc06_tp():
    # ecc03
    x = np.array([])
    y = np.array([])
    x = np.append(x, handle1.tot_file_shift['ecc06_1_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_1_y1y'])                                # Plot detail. DO NOT delete.!!!!
    x = np.append(x, handle1.tot_file_shift['ecc06_2_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_2_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc06_3_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_3_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc06_4_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_4_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc06_5_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_5_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc06_6_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_6_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc06_7_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc06_7_y1y'])
    # x = x - np.mean(x)
    print('tot_pos_overlay_shift ecc06:' + str(np.mean(x)))
    return x, y
def ecc08_tp():
    x = np.array([])
    y = np.array([])
    x = np.append(x, handle1.tot_file_shift['ecc08_1_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_1_y1y'])  # Plot detail. DO NOT delete.!!!!
    x = np.append(x, handle1.tot_file_shift['ecc08_2_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_2_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc08_3_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_3_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc08_4_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_4_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc08_5_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_5_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc08_6_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_6_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc08_7_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_7_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc08_8_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc08_8_y1y'])
    # x = x - np.mean(x)
    print('tot_pos_overlay_shift ecc08:' + str(np.mean(x)))
    return x, y
def ecc09_tp():
    x = np.array([])
    y = np.array([])
    x = np.append(x, handle1.tot_file_shift['ecc09_1_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_1_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_2_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_2_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_3_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_3_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_4_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_4_y1y'])                             # Plot detail. DO NOT delete.!!!!
    x = np.append(x, handle1.tot_file_shift['ecc09_5_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_5_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_6_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_6_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_7_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_7_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_8_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_8_y1y'])
    x = np.append(x, handle1.tot_file_shift['ecc09_9_y1x'])
    y = np.append(y, handle1.tot_file_shift['ecc09_9_y1y'])
    return x, y
def ecc095_tp():
    x1 = np.array([])
    y1 = np.array([])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_1_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_1_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_2_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_2_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_3_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_3_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_4_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_4_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_5_y1x'][:3000])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_5_y1y'][:3000])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_6_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_6_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_7_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_7_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_8_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_8_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_9_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_9_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_10_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_10_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_11_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_11_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_12_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_12_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_13_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_13_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_14_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_14_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc095_15_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc095_15_y1y'])
    return  x1, y1
def ecc098_tp():
    x1 = np.array([])
    y1 = np.array([])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_1_y1x']) #1
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_1_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_2_y1x']) #1
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_2_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_3_y1x']) #1
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_3_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_4_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_4_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_5_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_5_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_6_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_6_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_8_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_8_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc098_9_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc098_9_y1y'])
    # x1 = x1 - np.mean(x1)
    print('tot_pos_overlay_shift ecc098:'+str(np.mean(x1)))
    return x1, y1
def ecc0995_tp():
# ecc0995(T4-plasmid)
    x1 = np.array([])
    y1 = np.array([])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_1_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_1_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_2_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_2_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_3_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_3_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_4_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_4_y1y'])
    # x1 = np.append(x1, handle1.tot_file_shift['ecc0995_5_y1x']) # Cavity can't be aligned with T4 location due to T4 bleeching. Can't be included.
    # y1 = np.append(y1, handle1.tot_file_shift['ecc0995_5_y1y'])
    # x1 = np.append(x1, handle1.tot_file_shift['ecc0995_6_y1x'])
    # y1 = np.append(y1, handle1.tot_file_shift['ecc0995_6_y1y'])

    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_7_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_7_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_8_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_8_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_9_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_9_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_10_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_10_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_11_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_11_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_12_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_12_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_13_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_13_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_14_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_14_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_15_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_15_y1y'])
    x1 = np.append(x1, handle1.tot_file_shift['ecc0995_16_y1x'])
    y1 = np.append(y1, handle1.tot_file_shift['ecc0995_16_y1y'])
    return x1, y1

## Cavity overlay
def cav_ecc06():
    a, b = module.ellipse_para(ecc=0)
    n = 1000 # sampling pts
    t = np.linspace(0, 2*np.pi, n)
    x = a * np.cos(t)
    y = b * np.sin(t)
    return x, y

###Data setup##########
x0, y0 = ecc0_tp()
x1, y1 = ecc06_tp()
x2, y2 = ecc08_tp()
x3, y3 = ecc09_tp()
x4, y4 = ecc095_tp()
x5, y5 = ecc098_tp()
x6, y6 = ecc0995_tp()
X = [x0, x1, x2, x3, x4, x5, x6] # Input data array
Y = [y0, y1, y2, y3, y4, y5, y6]
# savename = 'Ecc0995_p.eps'
savename_tot = ['Ecc0_p.png', 'Ecc06_p.png', 'Ecc08_p.png', 'Ecc09_p.png',
            'Ecc095_p.png', 'Ecc098_p.png', 'Ecc0995_p.png']
figtitle_tot = ['Plasmid position distribution Ecc=0', 'Plasmid position distribution Ecc=0.6',
            'Plasmid position distribution Ecc=0.8', 'Plasmid position distribution Ecc=0.9',
            'Plasmid position distribution Ecc=0.95', 'Plasmid position distribution Ecc=0.98',
            'Plasmid position distribution Ecc=0.995']

savepath = '/media/zezhou/Se' \
           'agate Expansion Drive/McGillResearch/2019Manuscript_Analysis/Analysis/Plots/tplasmid/hist/'
xlim_tot = [10, 10, 10, 12, 12, 16, 16]
scalebar_tot = np.array(xlim_tot)/10*6.25 # Keep length of the scale bar the same

# Annotation text
marker_tot = ['e=0', 'e=0.6', 'e=0.8',
              'e=0.9', 'e=0.95', 'e=0.98', 'e=0.995']
scalelabel_tot = ['1um', '1um', '1um', '1.2um', '1.2um', '1.6um', '1.6um']

###
bins = 70
cmap = 'inferno'
fontsize = 25
labelsize = 12
cb_fontsize = 12
###############
a = 10
fig2 = plt.figure(figsize=[a, a/1.3333]) # figs8ize=[6.4, 4.8]

#### Do NOT change here. This section changes color ONLY ####
vmax = 0
vmin = 1000
############################
for i in range(7):
    x = X[i]
    y = Y[i]
    xlim = xlim_tot[i]

    # # Rescaling
    # h = np.histogram2d(x, y, bins=[bins, bins], range=[[-xlim, xlim], [-xlim, xlim]], density=True)
    # ## scale the color
    # sh = h[0][int(bins/2-10):int(bins/2+10), int(bins/2-4):int(bins/2+4)]
    # mi = np.min(sh) # smallest probability of the center square
    # ma = np.max(h[0]) #largest
    # ####################
    # if vmin>mi:
    #     vmin = mi
    # if vmax<ma:
    #     vmax = ma

width = 0.25
os.chdir("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/")
sfname = ['Ecc0', 'Ecc06', 'Ecc08',
              'Ecc09', 'Ecc095', 'Ecc098', 'Ecc0995']
vmax = 0
for i in range(7):
    x = X[i]
    y = Y[i]
    xlim = xlim_tot[i]
    figtitle = figtitle_tot[i]
    savename = savename_tot[i]
    # scale = 6.4/4.8
    c = i//4
    r = i%4
    ax2 = fig2.add_axes([0+c*width,1-(r+1)*width,width,width]) # Axes location. Right now it's an easy going version.
    # 2d hist

    # individual plot: uncomment this section for individual colorbar plot
    # h, xedges, yedges, img = ax2.hist2d(x, y, bins=[bins, bins], range=[[-xlim, xlim], [-xlim, xlim]], cmap=cmap, vmin=np.min(sh),
    #                                     vmax=np.max(h[0]), density=True, norm=mcolors.PowerNorm(0.7))
    weight = np.ones_like(x)/len(x)
    # Universal plot:
    if i !=0:
        h, xedges, yedges, img = ax2.hist2d(x, y, bins=[bins, bins], range=[[-xlim, xlim], [-xlim, xlim]], cmap=cmap,
                                            weights=weight)
    else:
        h, xedges, yedges, img = ax2.hist2d(x, y, bins=[bins, bins], range=[[-xlim, xlim], [-xlim, xlim]], cmap=cmap,
                                            norm = mcolors.PowerNorm(0.6) ,weights=weight)
        np.savetxt("Ecc0x_add.txt", xedges)
        np.savetxt("Ecc0y_add.txt", yedges)
        np.savetxt("Ecc0h_add.txt", h)
    print('current vmax:'+str(np.max(h)))
    print('current norm:'+str(np.sum(h)))
    if i==0:
        vmax=np.max(h)
    # colorbar: uncomment this section for individual colorbar
    # norm = colors.Normalize(vmin=np.min(h), vmax=np.max(h))
    # cb = fig2.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax = ax2)
    # cb.set_label('Probability', fontsize = fontsize, rotation = -90, horizontalalignment = 'center', verticalalignment = 'bottom')
    # sfntmp = sfname[i]

    # axes label
    ax2.set_xlabel(r'Position($\mu$m)', fontsize = fontsize)
    ax2.set_ylabel(r'Position($\mu$m)', fontsize = fontsize)

    # tick location and tick label
    xtick_temp = np.arange(0, xlim, step = 6.25)
    xtick_temp = np.concatenate((-xtick_temp, xtick_temp))
    xtick_temp = np.sort(xtick_temp)
    xtick_temp = np.delete(xtick_temp, len(xtick_temp)/2 - 1 )
    xticks = list(xtick_temp)
    xticklabel = []
    for ticks in xticks:
        xticklabel.append(str(ticks/6.25))
    ax2.set_xticks(xticks)
    ax2.set_yticks(xticks)

    # Tick label. Uncomment below to show unit in um
    ax2.set_xticklabels(xticklabel, fontsize = labelsize)
    ax2.set_yticklabels(xticklabel, fontsize = labelsize)

    # Axis switch
    ax2.axis('off')
    # Grid
    # ax2.grid(b=True, ls = '--', dash_capstyle='round')

    # Scale bar

    # Cavity overlay
    # ax2.plot(x, y, '-')

    # Title
    # ax2.set_title(figtitle, fontsize = fontsize)

    # text comment
    mark = marker_tot[i]
    ax2.text(x = -xlim*0.9, y=xlim*0.3, s=mark, fontsize = labelsize+7, color='w', fontweight='bold')

    # scale bar location
    sc = scalebar_tot[i]
    # rect = patches.Rectangle((xlim*0.3, -xlim*0.65), width = sc, height=sc/8, fill=True, color = 'white')
    rect = patches.Rectangle((0, -sc), width=sc, height=sc / 8, fill=True, color='white')
    ax2.add_patch(rect)
    # scale bar label
    # scalelabel = scalelabel_tot[i]
    # ax2.text(x=xlim * 0.3, y=-xlim * 0.95, s=scalelabel, fontsize=labelsize, color = 'white')
# fig2.suptitle('Plasmid distribution in different cavities', fontsize=20)

# Universal colorbar
cax = fig2.add_axes([1.2*width, 0.05*width, 0.03, 0.9*width])
norm = colors.Normalize(vmin=0, vmax=vmax)
print(vmax)
cb = fig2.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, orientation = 'vertical')

cb = fig2.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, orientation = 'vertical', ticks=[0,0.001,0.002])
# , ticks=[0,0.0025, 0.005]
cb.ax.tick_params(labelsize=cb_fontsize)
cb.ax.set_yticklabels(['0',r'$1e^{-3}$',  r'$2e^{-3}$'], fontdict={'fontsize': cb_fontsize+5})
cb.set_label('Prob.', fontsize = cb_fontsize+10 , horizontalalignment = 'center', rotation=-90, labelpad=20)
plt.show()

def crossection():
    fig = plt.figure(tight_layout=True)
    os.chdir("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/longaxes")
    areasq= []
    xdiv = []
    ydiv = []
    for i in range(7):
        ax = fig.add_subplot(4,2,i+1)
        x = X[i]
        y = Y[i]
        xlim = xlim_tot[i]
        weight = np.ones_like(x) / len(x)
        h, xedges, yedges= np.histogram2d(x, y, bins=[bins, bins], range=[[-xlim, xlim], [-xlim, xlim]],weights=weight)
        ctemp = int(bins/2)
        # short axes
        # ax.errorbar((xedges[:-1]+xedges[1:])/2, np.mean(h[(ctemp-2):(ctemp+1),: ],axis=0), yerr=np.std(h[:, (ctemp-2):(ctemp+1)],axis=1)/np.sqrt(3))
        # np.savetxt('ecc_'+str(i+1)+'_xedge.txt', (xedges[:-1]+xedges[1:])/2)
        # np.savetxt('ecc_' + str(i+1) + '_mean.txt', np.mean(h[(ctemp-2):(ctemp+1), :],axis=0))
        # np.savetxt('ecc_' + str(i+1) + '_yerr.txt', np.std(h[(ctemp-2):(ctemp+1), :],axis=0)/np.sqrt(3))

        # long axes
        ax.errorbar((xedges[:-1] + xedges[1:]) / 2, np.mean(h[:, (ctemp - 2):(ctemp + 1)], axis=1),
                    yerr=np.std(h[:, (ctemp - 2):(ctemp + 1)], axis=1) / np.sqrt(3))
        np.savetxt('ecc_' + str(i + 1) + '_xedge.txt', (xedges[:-1] + xedges[1:]) / 2)
        np.savetxt('ecc_' + str(i + 1) + '_mean.txt', np.mean(h[:, (ctemp - 2):(ctemp + 1)], axis=1))
        np.savetxt('ecc_' + str(i + 1) + '_yerr.txt', np.std(h[:, (ctemp - 2):(ctemp + 1)], axis=1) / np.sqrt(3))
        areasq.append((yedges[1]-yedges[0])*(xedges[1]-xedges[0]))
        xdiv.append((xedges[1]-xedges[0]))
        ydiv.append((yedges[1]-yedges[0]))
    np.savetxt('areasq.txt', areasq)
    np.savetxt('xdiv.txt', xdiv)
    np.savetxt('ydiv.txt', ydiv)
    # plt.savefig('histall_uni.png')
    plt.show()
    return
crossection()