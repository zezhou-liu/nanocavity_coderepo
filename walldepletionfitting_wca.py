from scipy.optimize import root
import os
import time
import module
import numpy as np
import matplotlib.pyplot as plt
import module
from scipy.optimize import minimize
from mpl_toolkits.mplot3d import Axes3D
import time
import matplotlib.colors as colors
import matplotlib.cm as cm
import scipy.integrate as integrate
import matplotlib.tri as tri
def distance_generation():
    filepath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/"
    os.chdir(filepath)

    # Change this to loop to process all data
    tmp = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/sampling/ecc0995_sampling.txt"
    data = np.loadtxt(tmp, skiprows=2)  # 1st column: x axis; 2nd column: y axis. 3rd column: z axis
    # Read the coordinate lattice
    x = data[:, 0]
    y = data[:, 1]
    a, b = module.ellipse_para(ecc=0.995)
    wall = np.zeros_like(x)
    # Create the depletion z profile
    for i in range(len(x)):
        xtmp = x[i]
        ytmp = y[i]
        rdis = xtmp ** 2 / a ** 2 + ytmp ** 2 / b ** 2
        if rdis > 1:
            continue

        def depletion(r):
            # r: dist to boundary; x,y: current coordinate of the point of interest; a,b: ecc's parameter
            return xtmp ** 2*(b - r) ** 2 + ytmp ** 2 * (a - r) ** 2 - (a - r) ** 2*(b - r) ** 2

        sol = root(depletion, x0=1e-5, method='lm')
        wall[i] = sol.x
    np.savetxt('ecc0995_rdist.txt', wall)
    plt.plot(wall, 'r')
    plt.show()
    return
def wdf(x):
    ## Using the following local variables to avoid loading everytime ## ADD THESE TO THE MAIN CODE BLOCK BEFORE USING THIS FUNCTION #####
    # rdistpath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/rdist/ecc0_rdist.txt"
    # simpath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/sampling/ecc0_sampling.txt"
    # simdata = np.loadtxt(simpath, skiprows=2)
    # a, b = module.ellipse_para(ecc=0)
    #########################################################################################################
    #######################
    # Initial guess value #
    A = x[0]
    rc = x[1]
    B = np.abs(x[2])
    #######################
    tgv=1e30
    # Set the outside probability equal to 0. Here I set t4 concentration to 1e30(tgv) to make sure the particle does not penetrate
    simdata[:, 2] = simdata[:, 2]/np.sum(simdata[:, 2])
    for i in range(len(simdata[:,0])):
        xtmp = simdata[i,0]
        ytmp = simdata[i,1]
        rdis = xtmp**2/a**2+ytmp**2/b**2
        if rdis >= 1.0:
            simdata[i,2]=tgv

    rdistdata = np.loadtxt(rdistpath)
    x_sim = np.array(simdata[:, 0]*100-simdata[0, 0]*100+0.04, dtype=int)
    y_sim = np.array(simdata[:, 1]*100-simdata[0, 1]*100-0.01, dtype=int)
    z_mesh = np.zeros([nx+1,ny+1])
    z_mesh[x_sim, y_sim] = A*np.exp(-rdistdata/np.abs(rc)) + abs(B)*simdata[:,2]
    z_mesh[-1,:] = tgv
    z_mesh[:,-1] = tgv
    prob = np.exp(-z_mesh)
    prob = np.transpose(prob/np.sum(prob))
    # plt.imshow(prob)
    # plt.show()
    return prob
def simprob_integrationver(x):
    A = x[0]
    rdc = x[1]
    B = np.abs(x[2])
    #######################
    tgv = 1e12
    # Set the outside probability equal to 0. Here I set t4 concentration to 1e30(tgv) to make sure the particle does not penetrate
    a,b = module.ellipse_para(ecc)
    outmask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2)>=1
    ured[outmask] = tgv
    wdp_intver = potential_generation(rdc) # wall depletion potential. Integration version
    landscape = A*wdp_intver+B*ured
    np.savetxt('landscape_ecc0_add.txt', landscape)
    prob = np.exp(-landscape)
    prob = np.transpose(prob / np.sum(prob))
    # plt.tricontourf(xred[inmask], yred[inmask], prob[inmask], cmap='inferno')
    # plt.show()
    # plt.imshow(prob)
    # plt.show()
    return prob
def simprob_integrationver_wca(x):
    A = x[0]
    sig = x[1]
    B = np.abs(x[2])
    #######################
    tgv = 1e12
    # Set the outside probability equal to 0. Here I set t4 concentration to 1e30(tgv) to make sure the particle does not penetrate
    a,b = module.ellipse_para(ecc)
    outmask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2)>=1
    ured[outmask] = tgv
    wdp_intver = potential_generation_wca(sig) # wall depletion potential. Integration version
    landscape = A*wdp_intver+B*ured
    # os.chdir('/home/zezhou/Desktop/')
    # np.savetxt('landscape_ecc0.txt', landscape)
    prob = np.exp(-landscape)
    prob = np.transpose(prob / np.sum(prob))

    # plt.tricontourf(xred[inmask], yred[inmask], prob[inmask], cmap='inferno')
    # plt.show()
    # plt.imshow(prob)
    # plt.show()
    return prob
def feedthrough(prob, xs, ys, h, x, y):
    # Convert the size of simulation to real data size
    probout = np.zeros_like(h)
    # Edge to center
    xsm, ysm = np.meshgrid(xs, ys)

    for i in range(len(x)-1):
        xl = x[i]  # x left edge
        xr = x[i + 1]  # x right edge
        for j in range(len(y)-1):
            yb = y[j] # y bottom edge
            yt = y[j+1] # y top edge
            maskx = (xsm >= xl)*(xsm<xr)
            masky = (ysm >= yb) * (ysm < yt)
            mask = maskx * masky
            noft = max(np.sum(mask), 1)
            probout[j, i] = np.sum(prob[mask])/noft
    return probout/np.sum(probout)
def feedthrough_intver(prob, xs, ys, h, x, y):
    # Convert the size of simulation to real data size
    probout = np.zeros_like(h)
    for i in range(len(x)-1):
        xl = x[i]  # x left edge
        xr = x[i + 1]  # x right edge
        for j in range(len(y)-1):
            yb = y[j] # y bottom edge
            yt = y[j+1] # y top edge
            maskx = (xs >= xl)*(xs<xr)
            masky = (ys >= yb) * (ys < yt)
            mask = maskx * masky
            noft = max(np.sum(mask), 1)
            probout[j, i] = np.sum(prob[mask])/noft
    return probout/np.sum(probout)
def potentialintegration(x, y, e, rdc):
    # x,y: coordinate for the point of interest
    # r, theta: polar coordinate for the point of interest
    # e: eccentricity of the cavity
    # rdc: decay length for the wall depletion
    a, b = module.ellipse_para(e)
    def intfunc(theta):
        dx = a*np.cos(theta)-x
        dy = b*np.sin(theta)-y
        return np.exp(-np.sqrt(dx**2+dy**2)/rdc)*np.sqrt(a**2*(np.sin(theta))**2+b**2*(np.cos(theta))**2)
    # t = time.time()
    result = integrate.quad(intfunc,0,2*np.pi)
    # print('single int time:'+str(time.time()-t))
    return result[0]
def potentialintegration_wca(x, y, e, sig):
    # x,y: coordinate for the point of interest
    # r, theta: polar coordinate for the point of interest
    # e: eccentricity of the cavity
    # rdc: decay length for the wall depletion
    a, b = module.ellipse_para(e)
    tgv = 1e12
    def intfunc(theta):
        dx = a*np.cos(theta)-x
        dy = b*np.sin(theta)-y
        r = np.sqrt(dx**2+dy**2)
        if r < sig*2**(1/6):
            out = 4*((sig/r)**12-(sig/r)**6+1/4)*np.sqrt(a**2*(np.sin(theta))**2+b**2*(np.cos(theta))**2)
        else:
            out = 0

        if out>=tgv:
            out=tgv

        return out
    # t = time.time()
    result = integrate.quad(intfunc,0,2*np.pi,points=np.linspace(0,2*np.pi,10))
    # print('single int time:'+str(time.time()-t))
    return result[0]
def potential_generation_testversion():
    filepath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/"
    os.chdir(filepath)

    # Change this to loop to process all data
    tmp = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/sampling/ecc08_sampling.txt"
    data = np.loadtxt(tmp, skiprows=2)  # 1st column: x axis; 2nd column: y axis. 3rd column: z axis
    # Read the coordinate lattice
    e = 0.8
    rdc = 0.6
    x = data[:, 0]
    y = data[:, 1]
    reduce_samplemask = (np.arange(0, len(x))%samplescale==0) # using local variable here.
    xred = x[reduce_samplemask]
    yred = y[reduce_samplemask]
    a, b = module.ellipse_para(e)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2)<1
    z = []
    tgv = 1
    tmp = time.time()
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    print("Start integration calculating")
    print("Current time:"+str(tmp))
    for i in range(len(xred)):
        # if i%200==100:
        #     # print("Step:"+str(i)+" in "+ str(len(x)))
        #     # print("Time taking:"+str(time.time()-tmp))
        #     amin = np.min(z)
        #     amax= np.max(z)
        #     level = np.linspace(amin,amax,20).tolist()
        #     # tcf = ax1.tricontourf(xred[:i], yred[:i], z, cmap='inferno',levels=level)
        #     # plt.pause(0.01)
        if mask[i]==0:
            z.append(tgv)
        else:
            z.append(potentialintegration(xred[i], yred[i],e,rdc))

    # os.chdir("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/walldepletion/")
    # np.savetxt('wdp_'+str(e)+'_'+str(rdc)+'_'+str(samplescale)+'samplescale', np.array(z))
    print("Finish calculation, ellapsed time:" + str(time.time()-tmp))
    # tcf = ax1.tricontourf(xred[mask], yred[mask], np.array(z)[mask], cmap='inferno',levels=level)
    # ax1.set_xticks([-1, 0, 1])
    # ax1.set_yticks([-0.8, 0, 0.8])
    # plt.tight_layout()
    # plt.show()
    return np.array(z)
def potential_generation(rdc):
    # Read the coordinate lattice
    a, b = module.ellipse_para(ecc)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1
    z = []
    tgv = 1
    tmp = time.time()
    print("Start integration calculating")
    print("Current time:" + str(tmp))
    for i in range(len(xred)):
        if mask[i] == 0:
            z.append(tgv)
        else:
            z.append(potentialintegration(xred[i], yred[i], ecc, rdc))
    print("Finish calculation, ellapsed time:" + str(time.time() - tmp))
    return np.array(z)
def potential_generation_wca(sig):
    # Read the coordinate lattice
    a, b = module.ellipse_para(ecc)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1
    z = []
    tgv = 1
    tmp = time.time()
    print("Start integration calculating")
    print("Current time:" + str(tmp))
    for i in range(len(xred)):
        if mask[i] == 0:
            z.append(tgv)
        else:
            z.append(potentialintegration_wca(xred[i], yred[i], ecc, sig))
    print("Finish calculation, ellapsed time:" + str(time.time() - tmp))
    # triang = tri.Triangulation(xred, yred)
    # plt.tricontourf(triang, z)
    # plt.show()
    return np.array(z)
def potential_generation_plotuseonly(rdc):
    # Read the coordinate lattice
    a, b = module.ellipse_para(ecc)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1
    z = []
    tgv = 1
    tmp = time.time()
    print("Start integration calculating")
    print("Current time:" + str(tmp))
    for i in range(len(xred)):
        if mask[i] == 0:
            z.append(tgv)
        else:
            z.append(potentialintegration(xred[i], yred[i], ecc, rdc))
    print("Finish calculation, ellapsed time:" + str(time.time() - tmp))
    return xred, yred, ured, np.array(z)
def potential_generation_plotuseonly_wca(rdc):
    # Read the coordinate lattice
    a, b = module.ellipse_para(ecc)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1
    z = []
    tgv = 1
    tmp = time.time()
    print("Start integration calculating")
    print("Current time:" + str(tmp))
    for i in range(len(xred)):
        if mask[i] == 0:
            z.append(tgv)
        else:
            z.append(potentialintegration_wca(xred[i], yred[i], ecc, rdc))
    print("Finish calculation, ellapsed time:" + str(time.time() - tmp))
    return xred, yred, ured, np.array(z)
# Loss function definition
def lossfunc(ini):
    # Loss function define by 2-norm integration. Can be modified.
    prob = wdf(ini) # Input is initial guess. Sim data loading is manually loaded by changing the filepath in wdf.
    prob = feedthrough(prob, xs, ys, h, x, y)
    global wdf_pltflag
    # Check plot of discrete sim result and data
    if wdf_pltflag == 1:
        plt.imshow(prob/np.sum(prob), cmap='inferno', aspect='equal')
        print('simu prob sum:'+str(np.sum(prob)))
        plt.show()
        plt.imshow(h, cmap='inferno', aspect='equal')
        print('real data prob sum:' + str(np.sum(h)))
        plt.show()
        wdf_pltflag = 0
    prob_flat = prob.flatten()
    h_flat = h.flatten()
    err = -np.inner(prob_flat, h_flat)/np.sqrt(np.dot(prob_flat, prob_flat))/np.sqrt(np.dot(h_flat,h_flat))
    print('current err:'+str(err))
    print('current norm:'+str(np.sum(prob)))
    return err
def lossfunc_integrationver(ini):
    # Loss function define by 2-norm integration. Can be modified.
    prob = simprob_integrationver(ini)  # Input is initial guess. Sim data loading is manually loaded by changing the filepath in wdf.
    prob = feedthrough_intver(prob, xred, yred, h, x, y)
    global wdf_pltflag

    # Check plot of discrete sim result and data
    if wdf_pltflag == 1:
        plt.imshow(prob, cmap='inferno')
        print('simu prob sum:' + str(np.sum(prob)))
        plt.show()
        plt.imshow(h, cmap='inferno', aspect='equal')
        print('real data prob sum:' + str(np.sum(h)))
        plt.show()
        wdf_pltflag = 0
    prob_flat = prob.flatten()
    h_flat = h.flatten()
    err = -np.inner(prob_flat, h_flat) / np.sqrt(np.dot(prob_flat, prob_flat)) / np.sqrt(np.dot(h_flat, h_flat))
    global iterindex
    plt.scatter(iterindex, err)
    iterindex += 1
    plt.pause(0.00001)
    with open('iterationindex.txt', 'a') as f:
        f.write(str(iterindex)+'\n')
    f.close()

    with open('sol.txt', 'a') as f:
        f.writelines(str(ini)+'\n')
    f.close()

    with open('err.txt', 'a') as f:
        f.writelines(str(err)+'\n')
    f.close()
    return err
def lossfunc_integrationver_wca(ini):
    # Loss function define by 2-norm integration. Can be modified.
    prob = simprob_integrationver_wca(ini)  # Input is initial guess. Sim data loading is manually loaded by changing the filepath in wdf.
    prob = feedthrough_intver(prob, xred, yred, h, x, y)
    global wdf_pltflag

    # Check plot of discrete sim result and data
    if wdf_pltflag == 1:
        # triang = tri.Triangulation(xred, yred)
        # plt.tricontourf(triang, probf)

        plt.imshow(prob, cmap='inferno')
        print('simu prob sum:' + str(np.sum(prob)))
        plt.show()
        plt.imshow(h, cmap='inferno', aspect='equal')
        print('real data prob sum:' + str(np.sum(h)))
        plt.show()
        wdf_pltflag = 0
    prob_flat = prob.flatten()
    h_flat = h.flatten()
    err = -np.inner(prob_flat, h_flat) / np.sqrt(np.dot(prob_flat, prob_flat)) / np.sqrt(np.dot(h_flat, h_flat)) #cosine similarity
    # err = 10**4*np.sum((prob_flat-h_flat)**2)#least square displacement
    print('Current sol:'+str(ini))
    global iterindex
    plt.scatter(iterindex, err)
    iterindex += 1
    plt.pause(0.00001)
    with open('iterationindex_add.txt', 'a') as f:
        f.write(str(iterindex)+'\n')
    f.close()

    with open('sol_add.txt', 'a') as f:
        f.writelines(str(ini)+'\n')
    f.close()

    with open('err_add.txt', 'a') as f:
        f.writelines(str(err)+'\n')
    f.close()
    return err
# Fitting by scipy.optimize.minimize
def fitting(sol_ini):
    # Initial guess of the solution
    method = 'Nelder-Mead'
    tic = time.clock()
    sol = minimize(lossfunc, x0=sol_ini, method=method)
    toc = time.clock()
    print(sol)
    print(method+"Optimizing time:"+str(toc-tic))
    return
def fitting_intver(sol_ini):
    # Initial guess of the solution
    method = 'Nelder-Mead'
    tic = time.clock()
    sol = minimize(lossfunc_integrationver, x0=sol_ini, method=method)
    toc = time.clock()
    print(sol)
    print(method+"Optimizing time:"+str(toc-tic))
    return
def fitting_intver_wca(sol_ini):
    # Initial guess of the solution
    method = 'Nelder-Mead'
    tic = time.clock()
    sol = minimize(lossfunc_integrationver_wca, x0=sol_ini, method=method)
    toc = time.clock()
    print(sol)
    print(method+"Optimizing time:"+str(toc-tic))
    return
def plt_t4():
    # plt T4 concentration
    t4labelsize=18
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    tcf = ax1.tricontourf(simdata[:, 0], simdata[:, 1], simdata[:, 2]/np.max(simdata[:, 2]), cmap='inferno')
    cb = fig1.colorbar(tcf,ticks=[0,0.5,1], shrink = 0.8)
    cb.ax.tick_params(labelsize=t4labelsize)
    cb.set_label('Normalized T4\nconcentration', fontsize = t4labelsize+3 , horizontalalignment = 'center', rotation=-90, labelpad=40)
    ax1.set_xticks([-1, 0, 1])
    ax1.set_xticklabels([-1, 0, 1], fontdict={'fontsize':t4labelsize})
    ax1.set_yticks([-0.8, 0, 0.8])
    ax1.set_yticklabels([-0.8, 0, 0.8], fontdict={'fontsize': t4labelsize})
    ax1.set_xlabel('X-position '
                   r'$(\mu m)$', fontsize=t4labelsize+3)
    ax1.set_ylabel('Y-position '
                   r'$(\mu m)$', fontsize=t4labelsize+3)
    plt.tight_layout()
    plt.show()
    return
def wall_deplt_intver(ini):
    # plt wall_dep concentration
    prob = simprob_integrationver(ini)
    a, b = module.ellipse_para(ecc)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1

    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    t4labelsize = 14
    # tcf = ax1.tricontourf(simdata[:, 0], simdata[:, 1], A*np.exp(-rdistdata/tau), cmap='inferno')
    tcf = ax1.tricontourf(xred[mask], yred[mask], prob[mask], cmap='inferno')
    cb = fig1.colorbar(tcf,ticks=[0,0.5,1], shrink = 0.8)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Normalized potential of \nwall-depletion', fontsize = 14 , horizontalalignment = 'center',
                 rotation=-90, labelpad=50)
    ax1.set_xticks([-1, 0, 1])
    ax1.set_xticklabels([-1, 0, 1], fontdict={'fontsize':t4labelsize})
    ax1.set_yticks([-0.8, 0, 0.8])
    ax1.set_yticklabels([-0.8, 0, 0.8], fontdict={'fontsize': t4labelsize})
    ax1.set_xlabel('X-position '
                   r'$(\mu m)$', fontsize=t4labelsize+3)
    ax1.set_ylabel('Y-position '
                   r'$(\mu m)$', fontsize=t4labelsize+3)
    plt.tight_layout()
    plt.show()
    return xred, yred, prob
def wall_deplt_intver_wca(ini):
    # plt wall_dep concentration
    prob = simprob_integrationver_wca(ini)
    a, b = module.ellipse_para(ecc)
    mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1

    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    t4labelsize = 14
    # tcf = ax1.tricontourf(simdata[:, 0], simdata[:, 1], A*np.exp(-rdistdata/tau), cmap='inferno')
    tcf = ax1.tricontourf(xred[mask], yred[mask], prob[mask], cmap='inferno')
    cb = fig1.colorbar(tcf,ticks=[0,0.5,1], shrink = 0.8)
    cb.ax.tick_params(labelsize=14)
    cb.set_label('Normalized potential of \nwall-depletion', fontsize = 14 , horizontalalignment = 'center',
                 rotation=-90, labelpad=50)
    ax1.set_xticks([-1, 0, 1])
    ax1.set_xticklabels([-1, 0, 1], fontdict={'fontsize':t4labelsize})
    ax1.set_yticks([-0.8, 0, 0.8])
    ax1.set_yticklabels([-0.8, 0, 0.8], fontdict={'fontsize': t4labelsize})
    ax1.set_xlabel('X-position '
                   r'$(\mu m)$', fontsize=t4labelsize+3)
    ax1.set_ylabel('Y-position '
                   r'$(\mu m)$', fontsize=t4labelsize+3)
    plt.tight_layout()
    plt.show()
    return xred, yred, prob
# def wall_deplt():
#     # plt wall_dep concentration
#     os.chdir("/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/walldepletion/")
#     r = np.loadtxt('wdp_0.8_0.6_80sample')
#     mask =(r!=1)
#     A = 1
#     tau=0.3
#     t4labelsize=18
#     rdistdata = np.loadtxt(rdistpath)
#
#     reduce_samplemask = (np.arange(0, len(simdata[:, 0])) % 12 == 0)
#     xred = simdata[:, 0][reduce_samplemask]
#     yred = simdata[:, 1][reduce_samplemask]
#     mask = (xred ** 2 / a ** 2 + yred ** 2 / b ** 2) < 1
#
#     fig1, ax1 = plt.subplots()
#     ax1.set_aspect('equal')
#     # tcf = ax1.tricontourf(simdata[:, 0], simdata[:, 1], A*np.exp(-rdistdata/tau), cmap='inferno')
#     tcf = ax1.tricontourf(xred[mask], yred[mask], r[mask], cmap='inferno')
#     cb = fig1.colorbar(tcf,ticks=[0,0.5,1], shrink = 0.8)
#     cb.ax.tick_params(labelsize=t4labelsize)
#     cb.set_label('Normalized potential of \nwall-depletion', fontsize = t4labelsize+3 , horizontalalignment = 'center',
#                  rotation=-90, labelpad=50)
#     ax1.set_xticks([-1, 0, 1])
#     ax1.set_xticklabels([-1, 0, 1], fontdict={'fontsize':t4labelsize})
#     ax1.set_yticks([-0.8, 0, 0.8])
#     ax1.set_yticklabels([-0.8, 0, 0.8], fontdict={'fontsize': t4labelsize})
#     ax1.set_xlabel('X-position '
#                    r'$(\mu m)$', fontsize=t4labelsize+3)
#     ax1.set_ylabel('Y-position '
#                    r'$(\mu m)$', fontsize=t4labelsize+3)
#     plt.tight_layout()
#     plt.show()
#     return
# def fitres_plt():
#     fig = plt.figure(figsize=[10, 10/1.3333])
#     width = 0.145
#     ecc_list = [0,0.6,0.8,0.9,0.95,0.98,0.995]
#     import matplotlib.patches as patches
#     for i in range(len(xslist)):
#         global simdata,rdistpath,nx,ny,xs,ys,h,x,y,a,b #many paras introduced in wdf. Don't remove it. Otherwise variable in main code will be referred.
#         c = i // 8
#         r = i % 8
#         ax = fig.add_axes([0 + c * width, 1 - (r + 1) * width, width, width])
#         nx = nxlist[i] # nofp in x dir.
#         ny = nylist[i]
#         xs = np.linspace(-xslist[i],xslist[i],nx+1)
#         ys = np.linspace(-yslist[i],yslist[i],ny+1)
#         h = np.loadtxt(datapath_list[i])
#         h = np.transpose(h / np.sum(h))
#         x = np.loadtxt(xpath_list[i]) / scale  # xedge from hist
#         y = np.loadtxt(ypath_list[i]) / scale
#         a, b = module.ellipse_para(ecc=ecc_list[i])
#         simdata = np.loadtxt(simpath_list[i], skiprows=2)
#         rdistpath = rdistpath_list[i]
#         prob = wdf(solution_set[i])  # Input is initial guess. Sim data loading is manually loaded by changing the filepath in wdf.
#         prob = feedthrough(prob, xs, ys, h, x, y)
#         print('current vmax:'+str(np.max(prob)))
#         ax.imshow(prob, cmap='inferno',)
#         ax.axis('off')
#     # Universal colorbar
#     cb_fontsize = 15
#     cax = fig.add_axes([0.15, 0.25, 0.015, 0.5])
#     norm = colors.Normalize(vmin=0, vmax=0.0023427399554377835)
#     # cb = fig2.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax = cax, orientation = 'vertical')
#     cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='inferno'), cax=cax, orientation='vertical',
#                        ticks=[0, 0.001, 0.002])
#
#     cb.ax.tick_params(labelsize=cb_fontsize)
#     cb.ax.set_yticklabels(['0', r'$1e^{-3}$', r'$2e^{-3}$'], fontdict={'fontsize': cb_fontsize + 5})
#     cb.set_label('Probability', fontsize=cb_fontsize + 10, horizontalalignment='center', rotation=-90, labelpad=20)
#     plt.show()
#     return
# def saveall():
#     os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/samplinginfo/')
#     with open("fit_solution.txt", 'w') as f:
#         for item in solution_set:
#             for idex, para in enumerate(item):
#                 if para != item[-1]:
#                     f.write(str(para)+"\t")
#                 else:
#                     f.write(str(para)+'\n')
#         f.close()
#     with open('xstot.txt', 'w') as f:
#         for item in xslist:
#             f.write(str(item)+'\t')
#         f.close()
#     with open('ystot.txt', 'w') as f:
#         for item in yslist:
#             f.write(str(item) + '\t')
#         f.close()
#     with open('nxtot.txt', 'w') as f:
#         for item in nxlist:
#             f.write(str(item) + '\t')
#         f.close()
#     with open('nytot.txt', 'w') as f:
#         for item in nylist:
#             f.write(str(item) + '\t')
#         f.close()

if __name__=='__main__':
    datapath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc0995h.txt"
    xpath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc0995x.txt"
    ypath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/data/Ecc0995y.txt"
    simpath = "/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/t4concentration/sampling/ecc0995_sampling.txt"
    ecc = 0.995 # Eccentricity
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/ecc0995')
    # xsample = 1.29099 # simulation xlim
    # ysample = 0.774597		 # simulation ylim
    # nx = 259 # nofp in x dir.
    # ny = 155 # nofp in y dir.
    sol_ini = [1, 0.1, 5e3] # Initial guess of the solution
    wdf_pltflag = 1 # if want to see the check plt of simulation
    scale = 6.25 # 6.25pixel/um
    samplescale = 1 # smallest resolution: 10nm. To improve the speed, samplescale decrease the resolution by
    iterindex = 0
    # <samplescale> fold. For example, if 10 is inserted, the simulation resolution will be 100nm, which is 10nm*samplescale.


    # Real data & mesh
    h = np.loadtxt(datapath)
    h = np.transpose(h/np.sum(h)) # Probability normalization
    x = np.loadtxt(xpath)/scale # xedge from hist
    y = np.loadtxt(ypath)/scale # yedge from

    # Simulation data. Integration version.
    simdata = np.loadtxt(simpath, skiprows=2)
    reduce_samplemask = (np.arange(0, len(simdata[:,0]))%samplescale==0) # Improve the speed
    xred = simdata[:,0][reduce_samplemask]
    yred = simdata[:,1][reduce_samplemask]
    ured = (simdata[:,2]**2/np.sum(simdata[:,2]**2))[reduce_samplemask] # concentration normalization

    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/simulationdata_wca/')
    prob = simprob_integrationver_wca([1.09470430e+00,1.58976501e-01,3.08520413e+03])
    # simprob_integrationver([ 7.15255093e+02 , 3.22745813e-02, -9.20572459e+03])
    np.savetxt('ecc0995x.txt',xred)
    np.savetxt('ecc0995y.txt',yred)
    np.savetxt('ecc0995_prob.txt', prob)
    # Fitting
#
# prob = simprob_integrationver_wca([1.84248642e-02, 2.64722009e-01, -1.83871839e+03])
# prob = feedthrough_intver(prob, xred, yred, h, x, y)
# plt.imshow(prob, cmap='inferno')
# plt.show()
# plt.imshow(h, cmap='inferno')
# plt.show()