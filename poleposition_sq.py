import walldepletionfitting as wf
import module
import os
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.colors as colors


def simu_generate():
    # Generate simulation data in folder /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/simulationdata/
    exp_path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()
    sol = []
    ecclib = []

    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc','')
            eccn = float(eccn[0]+'.'+eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath+folder)
            with open('sol.txt', 'r') as f:
                for temp_sol in f:
                    pass
            f.close()
            temp_sol = re.sub('\[', '', temp_sol)
            temp_sol = re.sub(']\n', '', temp_sol).split(' ')
            for num in temp_sol:
                if num != '':
                    sol_container.append(float(num))
            sol.append(sol_container)
    print(ecclib)
    sol = np.vstack(sol) # solution with the same order as ecclib
    print(np.mean(np.abs(sol),axis=0))
    print(np.std(np.abs(sol),axis=0)/np.sqrt(7))
    os.chdir(sim_comparepath)
    xred, yred, prob = wf.wall_deplt_intver(sol[2,:]) # Modify from here to return(line40-43). Also mody walldepletionfitting.py file line451 to line 453
    np.savetxt('ecc_7_xred.txt', xred)
    np.savetxt('ecc_7_yred.txt', yred)
    np.savetxt('ecc_7_prob.txt', prob)
    return
def wall_generate():
    # Generate simulation data in folder /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/simulationdata/
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/wallprofile/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()

    sol = []
    ecclib = []

    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc','')
            eccn = float(eccn[0]+'.'+eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath+folder)
            with open('sol.txt', 'r') as f:
                for temp_sol in f:
                    pass
            f.close()
            temp_sol = re.sub('\[', '', temp_sol)
            temp_sol = re.sub(']\n', '', temp_sol).split(' ')
            for num in temp_sol:
                if num != '':
                    sol_container.append(float(num))
            sol.append(sol_container)
    print(ecclib)
    sol = np.vstack(sol) # solution with the same order as ecclib
    os.chdir(sim_comparepath)
    x,y,u,prob = wf.potential_generation_plotuseonly(sol[2,1]) # Modify from here to return(line81-85). Also mody walldepletionfitting.py file line451 to line 453
    np.savetxt('ecc_7_xred.txt', x)
    np.savetxt('ecc_7_yred.txt', y)
    np.savetxt('ecc_7_concenred.txt', u)
    np.savetxt('ecc_7_wall.txt', prob)
    return
def get_exprange(numind):
    # numind: order of the file
    # xrange: column direction, upper and lower limit for the data range
    # yrange: row direction,...............
    os.chdir(exppath)
    xdiv = np.loadtxt('xdiv.txt')[numind]
    xedge = np.loadtxt('ecc_'+str(numind)+'_xedge.txt')
    yrange = [1.5*xdiv, -1.5*xdiv]
    xrange = [np.max(xedge), np.min(xedge)]
    return xrange, yrange
def plt_exp(ax, numind):
    # ax: axes name. numind:data index
    os.chdir(exppath+'longaxes')
    areasq = np.loadtxt('areasq.txt')[numind-1]
    xedge = np.loadtxt('ecc_'+str(numind)+'_xedge.txt')
    if numind == 3:
        xedge = xedge+0.1
    elif numind == 2:
        xedge = xedge + 0.4
    elif numind == 6:
        xedge =xedge+1
    ymean = np.loadtxt('ecc_' + str(numind) + '_mean.txt')/(areasq/6.25/6.25)
    yerr = np.loadtxt('ecc_' + str(numind) + '_yerr.txt')/(areasq/6.25/6.25)
    ax.errorbar((xedge)/6.25, ymean, yerr=yerr, ls = '', marker='o',ms=4, capsize=3., color='r')
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xedge', (xedge)/6.25)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/ymean', ymean)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/yerr', yerr)
    # print('a')
    # ax.plot((xedge) / 6.25, ymean, 'o')
    return ax
def plt_exp_short(ax, numind):
    # ax: axes name. numind:data index
    os.chdir(exppath)
    areasq = np.loadtxt('areasq.txt')[numind-1]
    xedge = np.loadtxt('ecc_'+str(numind)+'_xedge.txt')
    if numind == 2:
        xedge = xedge-0.3
    ymean = np.loadtxt('ecc_' + str(numind) + '_mean.txt')/(areasq/6.25/6.25)
    yerr = np.loadtxt('ecc_' + str(numind) + '_yerr.txt')/(areasq/6.25/6.25)
    ax.errorbar((xedge)/6.25, ymean, yerr=yerr, ls = '', marker='o',ms=4, capsize=3., color='r', label='Exp. data',zorder=0)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xedge', (xedge)/6.25)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/ymean', ymean)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/yerr', yerr)
    print('a')
    return ax
def longprojection_exp_simu_samegraph():
    global exppath
    simupath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'

    # loading simulation data
    os.chdir(simupath)
    filename = os.listdir()
    print(filename)
    order = []
    xredtot = []
    yredtot = []
    probtot = []

    for file in filename:
        if 'ecc' in file:
            name = file.split('_')
            if name[1] not in order:
                order.append(name[1])
                if name[1] == 'ecc0':
                    xredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'xredadd.txt'))
                    yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yredadd.txt'))
                    probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'probadd.txt'))
                else:
                    xredtot.append(np.loadtxt(name[0]+'_'+name[1]+'_'+'xred.txt'))
                    yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yred.txt'))
                    probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'prob.txt'))

    print(order)

    scale = 6.25

    ind = 0
    os.chdir(exppath + 'longaxes')
    colorint = 37
    colorshift = 0.8
    fig = plt.figure(figsize=[5,12],tight_layout=True)
    elabel = ['0','0.6','0.8','0.9','0.95','0.98','0.995']
    width = 0.7
    height = 0.12
    for index in order:
        ax = fig.add_axes([0.2, 0.95-int(index)*height,width, height])
        numind = int(index)
        tcmap = cm.get_cmap("inferno").colors
        xdiv = np.loadtxt('xdiv.txt')[numind-1]/scale
        ydiv = np.loadtxt('ydiv.txt')[numind-1]/scale
        xedge_t = np.loadtxt('ecc_' + str(numind) + '_xedge.txt')/scale
        xedge = np.append((xedge_t - ydiv/2.), xedge_t[-1]+ydiv/2.)
        simu_mean = []
        simu_sem = []
        area = (np.max(xredtot[ind])-np.min(xredtot[ind]))*(np.max(yredtot[ind])-np.min(yredtot[ind]))
        # print(area)
        probtot_den = probtot[ind]
        for i in range(len(xedge)-1):
            maskx = (xredtot[ind]>=xedge[i])*(xredtot[ind]<xedge[i+1])
            masky = (yredtot[ind]>=-1.5*xdiv)*(yredtot[ind]<1.5*xdiv)
            mask = maskx*masky
            npts = len(xredtot[ind])
            areasq = np.loadtxt('areasq.txt')[numind-1]
            if np.sum(mask)!=0:
                simu_mean.append(np.mean(probtot_den[mask])*areasq/6.25/6.25*npts/area/(areasq/6.25/6.25))
                simu_sem.append(np.std(probtot_den[mask])/np.sqrt(np.sum(mask))*areasq/6.25/6.25*npts/area/(areasq/6.25/6.25))
            else:
                simu_mean.append(0)
                simu_sem.append(0)
        os.chdir(exppath + 'longaxes')
        ax = plt_exp(ax, numind)
        ax.errorbar(xedge_t, simu_mean, yerr=simu_sem,label='e='+elabel[int(index)-1])
        if index =='7':
            ax.set_xlabel(r'Position $\left( \mu m\right)$', fontsize=15)
        else:
            ax.tick_params(axis='x',labelbottom=False)
        plt.legend(loc='upper right')
        ax.set_xlim([-2.7,2.7])
        ind+=1
    fig.text(0.04,0.5,r'Prob. density $\left(\mu m^{-2}\right)$', va='center', rotation='vertical', fontsize=15)
    plt.show()
    # y axis: probability density/ unit: prob/um^2
    return
def shortprojection_exp_simu_samegraph():
    global exppath
    simupath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/shortaxes'

    # loading simulation data
    os.chdir(simupath)
    filename = os.listdir()
    print(filename)
    order = []
    xredtot = []
    yredtot = []
    probtot = []

    for file in filename:
        if 'ecc' in file:
            name = file.split('_')
            if name[1] not in order:
                order.append(name[1])
                if name[1] == 'ecc0':
                    xredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'xredadd.txt'))
                    yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yredadd.txt'))
                    probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'probadd.txt'))
                else:
                    xredtot.append(np.loadtxt(name[0]+'_'+name[1]+'_'+'xred.txt'))
                    yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yred.txt'))
                    probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'prob.txt'))

    print(order)

    scale = 6.25
    os.chdir(exppath)
    ind = 0

    colorint = 37
    colorshift = 0.8
    fig = plt.figure(figsize=[5,12],tight_layout=True)
    elabel = ['0','0.6','0.8','0.9','0.95','0.98','0.995']
    width = 0.7
    height = 0.12
    for index in order:
        ax = fig.add_axes([0.2, 0.95-int(index)*height,width, height])
        numind = int(index)
        tcmap = cm.get_cmap("inferno").colors
        xdiv = np.loadtxt('xdiv.txt')[numind-1]/scale
        ydiv = np.loadtxt('ydiv.txt')[numind-1]/scale
        xedge_t = np.loadtxt('ecc_' + str(numind) + '_xedge.txt')/scale
        xedge = np.append((xedge_t - ydiv/2.), xedge_t[-1]+ydiv/2.)
        simu_mean = []
        simu_sem = []
        area = (np.max(xredtot[ind])-np.min(xredtot[ind]))*(np.max(yredtot[ind])-np.min(yredtot[ind]))
        # print(area)
        probtot_den = probtot[ind]
        for i in range(len(xedge)-1):
            maskx = (xredtot[ind]>=-1.5*xdiv)*(xredtot[ind]<1.5*xdiv)
            masky = (yredtot[ind]>=xedge[i])*(yredtot[ind]<xedge[i+1])
            mask = maskx*masky
            npts = len(xredtot[ind])
            areasq = np.loadtxt('areasq.txt')[numind-1]
            if np.sum(mask)!=0:
                simu_mean.append(np.mean(probtot_den[mask])*areasq/6.25/6.25*npts/area/(areasq/6.25/6.25))
                simu_sem.append(np.std(probtot_den[mask])/np.sqrt(np.sum(mask))*areasq/6.25/6.25*npts/area/(areasq/6.25/6.25))
            else:
                simu_mean.append(0)
                simu_sem.append(0)
        ax = plt_exp_short(ax, numind)
        ax.errorbar(xedge_t, simu_mean, yerr=simu_sem,label='e='+elabel[int(index)-1])
        if index =='7':
            ax.set_xlabel(r'Position $\left( \mu m\right)$', fontsize=15)
        else:
            ax.tick_params(axis='x',labelbottom=False)
        plt.legend(loc='upper right')
        ax.set_xlim([-2.7,2.7])
        ind+=1
    fig.text(0.04,0.5,r'Prob. density$\left(\mu m^{-2}\right)$', va='center', rotation='vertical', fontsize=15)
    plt.show()
    return
def longpotential_pocket():
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/wallprofile/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()

    sol = []
    ecclib = []
    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc','')
            eccn = float(eccn[0]+'.'+eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath+folder)
            with open('sol.txt', 'r') as f:
                for temp_sol in f:
                    pass
            f.close()
            temp_sol = re.sub('\[', '', temp_sol)
            temp_sol = re.sub(']\n', '', temp_sol).split(' ')
            for num in temp_sol:
                if num != '':
                    sol_container.append(float(num))
            sol.append(sol_container)
    sol = np.vstack(sol)
    orderecc = np.argsort(ecclib)
    sol = sol[orderecc,:]
    order = [1,2,3,4,5,6,7]
    fig = plt.figure(figsize=[5,12],tight_layout=True)
    elabel = ['0','0.6','0.8','0.9','0.95','0.98','0.995']
    width = 0.7
    height = 0.12
    stripwidth = 0.05 # unit:um
    os.chdir(sim_comparepath)
    for index in order:
        ax = fig.add_axes([0.2, 0.95-int(index)*height,width, height])

        xred = np.loadtxt('ecc_'+str(index)+'_xred.txt') # unit:um
        yred = np.loadtxt('ecc_' + str(index) + '_yred.txt')
        u0 = np.loadtxt('ecc_' + str(index) + '_concenred.txt')
        wall = np.loadtxt('ecc_' + str(index) + '_wall.txt')
        masky = (yred>=-stripwidth)*(yred<stripwidth)
        xred = xred[masky]
        orxred = np.argsort(xred)
        xred = xred[orxred]
        wall = np.abs(wall[masky] * sol[index - 1][0])[orxred]
        u0 = np.abs(u0[masky] * sol[index - 1][2])[orxred]

        x = []
        w = []
        u = []
        tmp = 'start'
        for i in xred:
            if i != tmp:
                tmp = i
                tmask = (xred == i)
                x.append(i)
                w.append(np.mean(wall[tmask]))
                u.append(np.mean(u0[tmask]))
        ax.plot(x, w, ls='dashdot',label = 'W')
        ax.plot(x, u,ls='dashed', label='E')
        ax.plot(x, np.array(w)+np.array(u), label='W+E')
        ax.set_ylim([0,2])
        ax.set_xlim([-3.2,3.2])
        ax.set_yticks([0,1])
        if index==7:
            ax.set_xlabel(r'Position $\left(\mu m\right)$', fontsize=15)
        elif index==1:
            ax.legend(loc='upper right', frameon=False)
    fig.text(0.04,0.5,r'Potential landscape $(k_BT)$', va='center', rotation='vertical', fontsize=15)
    plt.show()
def shortpotential_pocket():
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/wallprofile/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()

    sol = []
    ecclib = []
    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc', '')
            eccn = float(eccn[0] + '.' + eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath + folder)
            with open('sol.txt', 'r') as f:
                for temp_sol in f:
                    pass
            f.close()
            temp_sol = re.sub('\[', '', temp_sol)
            temp_sol = re.sub(']\n', '', temp_sol).split(' ')
            for num in temp_sol:
                if num != '':
                    sol_container.append(float(num))
            sol.append(sol_container)
    sol = np.vstack(sol)
    orderecc = np.argsort(ecclib)
    sol = sol[orderecc, :]
    order = [1, 2, 3, 4, 5, 6, 7]
    fig = plt.figure(figsize=[5, 12], tight_layout=True)
    elabel = ['0', '0.6', '0.8', '0.9', '0.95', '0.98', '0.995']
    width = 0.7
    height = 0.12
    stripwidth = 0.05  # unit:um
    os.chdir(sim_comparepath)
    for index in order:
        ax = fig.add_axes([0.2, 0.95 - int(index) * height, width, height])

        xred = np.loadtxt('ecc_' + str(index) + '_xred.txt')  # unit:um
        yred = np.loadtxt('ecc_' + str(index) + '_yred.txt')
        u0 = np.loadtxt('ecc_' + str(index) + '_concenred.txt')
        wall = np.loadtxt('ecc_' + str(index) + '_wall.txt')
        maskx = (xred >= -stripwidth) * (xred < stripwidth)
        yred = yred[maskx]
        oryred = np.argsort(yred)
        yred = yred[oryred]
        wall = np.abs(wall[maskx] * sol[index - 1][0])[oryred]
        u0 = np.abs(u0[maskx] * sol[index - 1][2])[oryred]

        x = []
        w = []
        u = []
        tmp = 'start'
        for i in yred:
            if i != tmp:
                tmp = i
                tmask = (yred == i)
                x.append(i)
                w.append(np.mean(wall[tmask]))
                u.append(np.mean(u0[tmask]))
        ax.plot(x, w, ls='dashdot', label='W')
        ax.plot(x, u, ls='dashed', label='E')
        ax.plot(x, (np.array(w) + np.array(u)), label='W+E')
        ax.set_ylim([0, 2])
        ax.set_xlim([-1,1])
        ax.tick_params(axis='x', labelbottom=False)
        ax.set_yticks([0, 1])
        if index == 7:
            ax.set_xlabel(r'Position $\left(\mu m\right)$', fontsize=15)
            ax.tick_params(axis='x', labelbottom=True)
        elif index == 1:
            ax.legend(loc='upper right', frameon=False)
    fig.text(0.04, 0.5, r'Potential landscape ($k_BT$)', va='center', rotation='vertical', fontsize=15)
    plt.show()
def simugraph():
    simupath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'

    # loading simulation data
    os.chdir(simupath)
    filename = os.listdir()
    print(filename)
    order = []
    xredtot = []
    yredtot = []
    probtot = []

    for file in filename:
        if 'ecc' in file:
            name = file.split('_')
            if name[1] not in order:
                order.append(name[1])
                xredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'xred.txt'))
                yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yred.txt'))
                probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'prob.txt'))

    print(order)

    scale = 6.25
    os.chdir(exppath)

    colorint = 37
    colorshift = 0.8
    a = 10
    fig = plt.figure(figsize=[a, a/1.3333])
    # elabel = ['0', '0.6','0.8', '0.9', '0.95', '0.98', '0.995'] full data
    elabel = ['0', '0.8', '0.9', '0.95', '0.98', '0.995']# Show purpose. Remove 0.6 for clarity
    width = 0.25
    height = 0.12
    # ecc = ['0','06','08','09','095','098','0995'] full data
    ecc = ['0', '08', '09', '095', '098', '0995']
    sct = [1,1,1,1.2,1.2,1.6,1.6]
    ind = 0

    order.remove('2') # Show purpose. Remove 0.6 for clarity

    for index in order:
        numinx = int(index)
        if numinx>=2: # Show purpose . Remove 0.6 for clarity
            numinx=numinx-1
        c = (numinx-1)// 3
        r = (numinx-1) % 3
        ax = fig.add_axes([0 + c * width, 1 - (r + 1) * width, width, width])  # Axes location. Right now it's an easy going version.

        xred = xredtot[ind]
        yred = yredtot[ind]
        prob = probtot[ind]

        h = np.loadtxt("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc"+ecc[numinx-1]+"h.txt")
        h = np.transpose(h / np.sum(h))  # Probability normalization
        x = np.loadtxt("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc"+ecc[numinx-1]+"x.txt") / scale  # xedge from hist
        y = np.loadtxt("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc"+ecc[numinx-1]+"y.txt") / scale  # yedge from
        X, Y = np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2)
        prob_simu = wf.feedthrough_intver(prob, xred, yred, h, x, y)

        ax.pcolormesh(X,Y,prob_simu,cmap='inferno')
        ax.tick_params(bottom=False,left=False,top=False,right=False,labelbottom=False,labelleft=False)
        # rect = patches.Rectangle((xlim*0.3, -xlim*0.65), width = sc, height=sc/8, fill=True, color = 'white')
        sc = sct[numinx-1]
        rect = patches.Rectangle((0, -sc), width=sc, height=sc / 8, fill=True, color='white')
        ax.add_patch(rect)
        mark = elabel[numinx-1]
        xlim = np.max(x)
        print('summation of the pixel value:'+str(np.sum(prob_simu)))
        if numinx==1:
            xxlim = np.max(prob_simu)
        ax.text(x=-xlim * 0.9, y=xlim * 0.3, s="e="+mark, fontsize=20, color='w', fontweight='bold')
        ind +=1
    # Universal colorbar
    # cax = fig.add_axes([1.2*width, 0.05*width, 0.03, 0.9*width])
    # cax = fig.add_axes([2.2 * width, 1.2 * width, 0.03, 1.6 * width])
    cax = fig.add_axes([2.1 * width, 1.2 * width, 0.03, 2.6 * width])
    norm = colors.Normalize(vmin=0, vmax=xxlim)
    print('color max value:'+str(xxlim))
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='inferno'), cax = cax, orientation = 'vertical',ticks=[0,0.0005, 0.001,0.0015])

    cb.ax.tick_params(labelsize=20)
    cb.ax.set_yticklabels(['0', r'$5e^{-4}$', r'$1e^{-3}$', r'$1.5e^{-3}$'], fontdict={'fontsize': 15})
    cb.set_label('Prob.', fontsize = 25 , horizontalalignment = 'center', rotation=-90, labelpad=20)
    plt.show()
    return
def plot_fig5_bcde():
    # use ecc=0.9 as an example.
    global exppath
    simupath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'

    # loading simulation data
    os.chdir(simupath)
    filename = os.listdir()
    print(filename)
    order = []
    xredtot = []
    yredtot = []
    probtot = []

    for file in filename:
        if 'ecc' in file:
            name = file.split('_')
            if name[1] not in order:
                order.append(name[1])
                xredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'xred.txt'))
                yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yred.txt'))
                probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'prob.txt'))

    print(order)

    scale = 6.25
    os.chdir(exppath)
    ind = -2

    colorint = 37
    colorshift = 0.8
    fig = plt.figure( tight_layout=True)
    width = 0.4
    height = 0.4

    # Experiment cross-section/simulation layover
    # long-axes
    ax = fig.add_subplot(2,2,1)
    numind = 4
    tcmap = cm.get_cmap("inferno").colors
    xdiv = np.loadtxt('xdiv.txt')[numind-1] / scale
    ydiv = np.loadtxt('ydiv.txt')[numind-1] / scale
    xedge_t = np.loadtxt('ecc_' + str(numind) + '_xedge.txt') / scale
    xedge = np.append((xedge_t - ydiv / 2.), xedge_t[-1] + ydiv / 2.)
    simu_mean = []
    simu_sem = []
    area = (np.max(xredtot[ind]) - np.min(xredtot[ind])) * (np.max(yredtot[ind]) - np.min(yredtot[ind]))
    # print(area)
    probtot_den = probtot[ind]
    for i in range(len(xedge) - 1):
        maskx = (xredtot[ind] >= xedge[i]) * (xredtot[ind] < xedge[i + 1])
        masky = (yredtot[ind] >= -1.5 * xdiv) * (yredtot[ind] < 1.5 * xdiv)
        mask = maskx * masky
        npts = len(xredtot[ind])
        areasq = np.loadtxt('areasq.txt')[numind - 1]
        if np.sum(mask) != 0:
            simu_mean.append(np.mean(probtot_den[mask]) * areasq / 6.25 / 6.25 * npts / area / (areasq / 6.25 / 6.25))
            simu_sem.append(np.std(probtot_den[mask]) / np.sqrt(np.sum(mask)))
        else:
            simu_mean.append(0)
            simu_sem.append(0)
    ax = plt_exp(ax, numind)
    ax.errorbar(xedge_t, simu_mean, yerr=simu_sem, color='k', ls='--')
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xedge_t', xedge_t)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simu_mean', simu_mean)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simu_sem', simu_sem)
    ax.set_xlabel(r'Position $\left( \mu m\right)$', fontsize=15)
    ax.set_ylabel('Prob. density\n'
                  r'$\left(\mu m^{-2}\right)$', fontsize=15)
    ax.set_xlim([-1.7, 1.7])
    ax.set_ylim([-0.05, 0.9])
    # short-axes

    exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/shortaxes'

    ax = fig.add_subplot(2,2,2)
    numind = 4
    tcmap = cm.get_cmap("inferno").colors
    xdiv = np.loadtxt('xdiv.txt')[numind-1] / scale
    ydiv = np.loadtxt('ydiv.txt')[numind-1] / scale
    xedge_t = np.loadtxt('ecc_' + str(numind) + '_xedge.txt') / scale
    xedge = np.append((xedge_t - ydiv / 2.), xedge_t[-1] + ydiv / 2.)
    simu_mean = []
    simu_sem = []
    area = (np.max(xredtot[ind]) - np.min(xredtot[ind])) * (np.max(yredtot[ind]) - np.min(yredtot[ind]))
    # print(area)
    probtot_den = probtot[ind]
    for i in range(len(xedge) - 1):
        maskx = (xredtot[ind] >= -1.5 * xdiv) * (xredtot[ind] < 1.5 * xdiv)
        masky = (yredtot[ind] >= xedge[i]) * (yredtot[ind] < xedge[i + 1])
        mask = maskx * masky
        npts = len(xredtot[ind])
        areasq = np.loadtxt('areasq.txt')[numind - 1]
        if np.sum(mask) != 0:
            simu_mean.append(np.mean(probtot_den[mask]) * areasq / 6.25 / 6.25 * npts / area / (areasq / 6.25 / 6.25))
            simu_sem.append(np.std(probtot_den[mask]) / np.sqrt(np.sum(mask)))
        else:
            simu_mean.append(0)
            simu_sem.append(0)
    ax = plt_exp_short(ax, numind)
    ax.plot(xedge_t, simu_mean, label='Sim. fitting', color='k', ls='--',zorder=1)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/xedge_t', xedge_t)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simu_mean', simu_mean)
    # np.savetxt('/home/zezhou/McGillResearch/2019Manuscript_Analysis/Submission/ArchiveData/simu_sem', simu_sem)
    ax.set_xlabel(r'Position $\left( \mu m\right)$', fontsize=15)
    # ax.set_ylabel('Probability density\n'
    #               r'$\left(\mu m^{-2}\right)$', fontsize=15)
    ax.set_xlim([-1.7, 1.7])
    ax.set_ylim([-0.05, 0.9])
    ax.legend(frameon=False)
    # formation of the pocket
    ax=fig.add_subplot(223)
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/wallprofile/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()

    sol = []
    ecclib = []
    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc','')
            eccn = float(eccn[0]+'.'+eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath+folder)
            with open('sol.txt', 'r') as f:
                for temp_sol in f:
                    pass
            f.close()
            temp_sol = re.sub('\[', '', temp_sol)
            temp_sol = re.sub(']\n', '', temp_sol).split(' ')
            for num in temp_sol:
                if num != '':
                    sol_container.append(float(num))
            sol.append(sol_container)
    sol = np.vstack(sol)
    orderecc = np.argsort(ecclib)
    sol = sol[orderecc,:]
    order = [1,2,3,4,5,6,7]

    stripwidth = 0.05 # unit:um
    os.chdir(sim_comparepath)
    index=4
    xred = np.loadtxt('ecc_'+str(index)+'_xred.txt') # unit:um
    yred = np.loadtxt('ecc_' + str(index) + '_yred.txt')
    u0 = np.loadtxt('ecc_' + str(index) + '_concenred.txt')
    wall = np.loadtxt('ecc_' + str(index) + '_wall.txt')
    masky = (yred>=-stripwidth)*(yred<stripwidth)
    xred = xred[masky]
    orxred = np.argsort(xred)
    xred = xred[orxred]
    wall = np.abs(wall[masky] * sol[index - 1][0])[orxred]
    u0 = np.abs(u0[masky] * sol[index - 1][2])[orxred]

    x = []
    w = []
    u = []
    tmp = 'start'
    for i in xred:
        if i != tmp:
            tmp = i
            tmask = (xred == i)
            x.append(i)
            w.append(np.mean(wall[tmask]))
            u.append(np.mean(u0[tmask]))
    ax.plot(x, w, ls='dashdot',label = 'Wall potential',color='r')
    ax.plot(x, u,ls='dashed', label='Exclusive potential', color='k')
    ax.plot(x, np.array(w)+np.array(u), label='W+E', color='g')
    ax.set_ylabel(r'$U_{p}$', fontsize=15)
    ax.set_xlim([-1.7, 1.7])
    ax.set_ylim([-0.05, 1.5])
    ax.set_xlabel(r'Position $\left( \mu m\right)$', fontsize=15)
    # ax.legend()

    # zoom in
    ax=fig.add_axes([0.185,0.3,0.15,0.15])
    ax.tick_params(labelleft=False, labelright=True,direction='in',right=True,left=False)
    ax.plot(x, w, ls='dashdot', color='r')
    ax.plot(x, u,ls='dashed', color='k')
    ax.plot(x, np.array(w)+np.array(u),color='g')
    ax.set_xlim([1,1.3])
    ax.set_xticks([1,1.3])
    ax.set_ylim([-0.02,0.25])
    ax.set_yticks([0,0.2])

    # short axes
    ax=fig.add_subplot(224)
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/wallprofile/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()

    sol = []
    ecclib = []
    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc', '')
            eccn = float(eccn[0] + '.' + eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath + folder)
            with open('sol.txt', 'r') as f:
                for temp_sol in f:
                    pass
            f.close()
            temp_sol = re.sub('\[', '', temp_sol)
            temp_sol = re.sub(']\n', '', temp_sol).split(' ')
            for num in temp_sol:
                if num != '':
                    sol_container.append(float(num))
            sol.append(sol_container)
    sol = np.vstack(sol)
    orderecc = np.argsort(ecclib)
    sol = sol[orderecc, :]
    order = [1, 2, 3, 4, 5, 6, 7]

    stripwidth = 0.05  # unit:um
    os.chdir(sim_comparepath)
    index=4
    xred = np.loadtxt('ecc_' + str(index) + '_xred.txt')  # unit:um
    yred = np.loadtxt('ecc_' + str(index) + '_yred.txt')
    u0 = np.loadtxt('ecc_' + str(index) + '_concenred.txt')
    wall = np.loadtxt('ecc_' + str(index) + '_wall.txt')
    maskx = (xred >= -stripwidth) * (xred < stripwidth)
    yred = yred[maskx]
    oryred = np.argsort(yred)
    yred = yred[oryred]
    wall = np.abs(wall[maskx] * sol[index - 1][0])[oryred]
    u0 = np.abs(u0[maskx] * sol[index - 1][2])[oryred]

    x = []
    w = []
    u = []

    for i in yred:
        if i != tmp:
            tmp = i
            tmask = (yred == i)
            x.append(i)
            w.append(np.mean(wall[tmask]))
            u.append(np.mean(u0[tmask]))
    ax.plot(x, w, ls='dashdot', label='W', color='r')
    ax.plot(x, u, ls='dashed', label='E', color='k')
    ax.plot(x, np.array(w) + np.array(u), label='W+E', color='g')
    ax.set_xlim([-1.7, 1.7])
    ax.set_ylim([-0.05, 3])
    ax.set_xlabel(r'Position $\left(\mu m\right)$', fontsize=15)
    ax.legend(loc=(0.67,0.55), frameon=False,handlelength=1.5)

    # zoom in
    ax=fig.add_axes([0.625,0.3,0.15,0.15])
    ax.tick_params(labelleft=False, labelright=True,direction='in',right=True,left=False)
    ax.plot(x, w, ls='dashdot', color='r')
    ax.plot(x, u, ls='dashed',  color='k')
    ax.plot(x, np.array(w) + np.array(u), color='g')
    ax.set_xlim([0.3, 0.5])
    ax.set_xticks([0.3,0.5])
    ax.set_ylim([-0.05, 0.7])
    ax.set_yticks([0, 0.5])
    plt.show()
    return
def simu_generate_ecc0add():
    # Generate simulation data in folder /home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/simulationdata/
    exp_path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'
    sim_resultpath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/'
    sim_comparepath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    os.chdir(sim_resultpath)
    foldername = os.listdir()
    sol = []
    ecclib = []

    for folder in foldername:
        if 'ecc' in folder:
            sol_container = []  # save solution temporarily
            eccn = folder.replace('ecc','')
            eccn = float(eccn[0]+'.'+eccn[1:])
            ecclib.append(eccn)
            # cd into folder and read the converged solution
            os.chdir(sim_resultpath+folder)
            if folder=='ecc0':
                with open('sol_add.txt', 'r') as f:
                    for temp_sol in f:
                        pass
                f.close()
                temp_sol = re.sub('\[', '', temp_sol)
                temp_sol = re.sub(']\n', '', temp_sol).split(' ')
                for num in temp_sol:
                    if num != '':
                        sol_container.append(float(num))
            else:
                with open('sol.txt', 'r') as f:
                    for temp_sol in f:
                        pass
                f.close()
                temp_sol = re.sub('\[', '', temp_sol)
                temp_sol = re.sub(']\n', '', temp_sol).split(' ')
                for num in temp_sol:
                    if num != '':
                        sol_container.append(float(num))
            sol.append(sol_container)
    print(ecclib)
    sol = np.vstack(sol) # solution with the same order as ecclib
    print(np.mean(np.abs(sol),axis=0))
    print(np.std(np.abs(sol),axis=0)/np.sqrt(7))
    os.chdir(sim_comparepath)
    xred, yred, prob = wf.wall_deplt_intver(sol[-1,:]) # Modify from here to return(line40-43). Also mody walldepletionfitting.py file line451 to line 453

    np.savetxt('ecc_1_xredadd.txt', xred)
    np.savetxt('ecc_1_yredadd.txt', yred)
    np.savetxt('ecc_1_probadd.txt', prob)
    return
def simugraph_ecc0add():
    simupath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata/'
    exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'

    # loading simulation data
    os.chdir(simupath)
    filename = os.listdir()
    print(filename)
    order = []
    xredtot = []
    yredtot = []
    probtot = []

    for file in filename:
        if 'ecc' in file:
            name = file.split('_')
            if name[1] not in order:
                if name[1] == '1':
                    order.append(name[1])
                    xredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'xredadd.txt'))
                    yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yredadd.txt'))
                    probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'probadd.txt'))
                else:
                    order.append(name[1])
                    xredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'xred.txt'))
                    yredtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'yred.txt'))
                    probtot.append(np.loadtxt(name[0] + '_' + name[1] + '_' + 'prob.txt'))

    print(order)

    scale = 6.25
    os.chdir(exppath)

    colorint = 37
    colorshift = 0.8
    a = 10
    fig = plt.figure(figsize=[a, a/1.3333])
    # elabel = ['0', '0.6','0.8', '0.9', '0.95', '0.98', '0.995'] full data
    elabel = ['0', '0.8', '0.9', '0.95', '0.98', '0.995']# Show purpose. Remove 0.6 for clarity
    width = 0.25
    height = 0.12
    # ecc = ['0','06','08','09','095','098','0995'] full data
    ecc = ['0', '08', '09', '095', '098', '0995']
    sct = [1,1,1,1.2,1.2,1.6,1.6]
    ind = 0

    order.remove('2') # Show purpose. Remove 0.6 for clarity

    for index in order:
        numinx = int(index)
        if numinx>=2: # Show purpose . Remove 0.6 for clarity
            numinx=numinx-1
        c = (numinx-1)// 3
        r = (numinx-1) % 3
        ax = fig.add_axes([0 + c * width, 1 - (r + 1) * width, width, width])  # Axes location. Right now it's an easy going version.

        xred = xredtot[ind]
        yred = yredtot[ind]
        prob = probtot[ind]
        if ecc[numinx-1]=='0':
            h = np.loadtxt(
                "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc" + ecc[
                    numinx - 1] + "h_add.txt")
            h = np.transpose(h / np.sum(h))  # Probability normalization
            x = np.loadtxt(
                "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc" + ecc[
                    numinx - 1] + "x_add.txt") / scale  # xedge from hist
            y = np.loadtxt(
                "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc" + ecc[
                    numinx - 1] + "y_add.txt") / scale  # yedge from
        else:
            h = np.loadtxt("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc"+ecc[numinx-1]+"h.txt")
            h = np.transpose(h / np.sum(h))  # Probability normalization
            x = np.loadtxt("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc"+ecc[numinx-1]+"x.txt") / scale  # xedge from hist
            y = np.loadtxt("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc"+ecc[numinx-1]+"y.txt") / scale  # yedge from
        X, Y = np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2)
        prob_simu = wf.feedthrough_intver(prob, xred, yred, h, x, y)

        ax.pcolormesh(X,Y,prob_simu,cmap='inferno')
        ax.tick_params(bottom=False,left=False,top=False,right=False,labelbottom=False,labelleft=False)
        # rect = patches.Rectangle((xlim*0.3, -xlim*0.65), width = sc, height=sc/8, fill=True, color = 'white')
        sc = sct[numinx-1]
        rect = patches.Rectangle((0, -sc), width=sc, height=sc / 8, fill=True, color='white')
        ax.add_patch(rect)
        mark = elabel[numinx-1]
        xlim = np.max(x)
        print('summation of the pixel value:'+str(np.sum(prob_simu)))
        if numinx==1:
            xxlim = np.max(prob_simu)
        ax.text(x=-xlim * 0.9, y=xlim * 0.3, s="e="+mark, fontsize=20, color='w', fontweight='bold')
        ind +=1
    # Universal colorbar
    # cax = fig.add_axes([1.2*width, 0.05*width, 0.03, 0.9*width])
    # cax = fig.add_axes([2.2 * width, 1.2 * width, 0.03, 1.6 * width])
    cax = fig.add_axes([2.1 * width, 1.2 * width, 0.03, 2.6 * width])
    norm = colors.Normalize(vmin=0, vmax=xxlim)
    print('color max value:'+str(xxlim))
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='inferno'), cax = cax, orientation = 'vertical',ticks=[0,0.0005, 0.001,0.0015])

    cb.ax.tick_params(labelsize=20)
    cb.ax.set_yticklabels(['0', r'$5e^{-4}$', r'$1e^{-3}$', r'$1.5e^{-3}$'], fontdict={'fontsize': 15})
    cb.set_label('Prob.', fontsize = 25 , horizontalalignment = 'center', rotation=-90, labelpad=20)
    plt.show()
    return

shortprojection_exp_simu_samegraph()
# exppath = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata/experimentaldata/'
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# plt_exp(ax,3)
# ax.set_xlabel(r'Position $\left( \mu m\right)$', fontsize=15)
# ax.set_ylabel(r'Prob. density $\left(\mu m^{-2}\right)$', fontsize=15)
# ax.text(x=1,y=0.6, s='ecc = 0.8', fontsize=12)
# plt.show()