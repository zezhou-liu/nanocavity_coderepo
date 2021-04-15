import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time
import os
################# Parameter setup ###############
n = 5000 # Step number
n_test = 250
n_teststep = 5
l = 1e-1 # Step length

os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/TOTPOTENTIAL/Exponential/ecc0/')
################# Load potential file ###########
xland = np.loadtxt('x_add.txt')
yland = np.loadtxt('y_add.txt')
zland = np.loadtxt('landscape_ecc0_add.txt')
################# Generate traj #################

def brownian(n):
    pos = [] # Position container
    r1 = np.random.rand()
    r2 = np.random.rand()
    tt = np.random.rand()*2*np.pi
    inipos = np.array([r1*np.cos(tt), r2*np.sin(tt)])  # Initial position
    pos.append(inipos)

    theta = 2*np.pi*np.random.rand(n_teststep*n) # direction of the particle moving
    jumpprob = np.random.rand(n) # random number to determine if the particle is going to move or not

    # Setup the initial potential
    # pot_tmp = interpolate.griddata((xland, yland), zland, inipos, method='nearest')

    for i in range(n):
        theta_temp = theta[i*n_teststep:(i+1)*n_teststep]
        jump_temp = jumpprob[i]

        dx = l*np.cos(theta_temp)
        dy = l*np.sin(theta_temp)
        newpos = pos[-1]+np.transpose(np.array([dx,dy])) # test steps

        pot_new = interpolate.griddata((xland, yland), zland, newpos, method='nearest') # find the potential of those test steps

        Z = np.sum(np.exp(-pot_new)) # calculate partition function
        prob_array = np.cumsum(np.exp(-pot_new))/Z

        pos_idx = np.argmax(jump_temp<=prob_array) # find the selected test step
        pos.append(newpos[pos_idx]) # append the position
    pos = np.array(pos) # pos[x,y]
    return pos
def run():
    pos = []
    t = time.time()
    print('Start @:'+str(t))
    os.chdir('/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/ecc0/')
    file = open('pos3.txt','a') # appending mode
    for i in range(n_test):
        ts = time.time()
        temp_pos = brownian(n)
        ta = time.time()
        print('Elapsed time:'+str(ta-t)+' s')
        print('Estimate time left'+str((n_test-i)*(ta-ts))+' s')
        np.savetxt(file, temp_pos) # appending new position traj to the same file.
def plotter():
    x = np.loadtxt('/home/zezhou/Desktop/pos4.txt')
    plt.hist2d(x[:,0], x[:,1], bins=20)
    plt.show()

run()