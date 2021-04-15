import numpy as np
import os

path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exp/MSDfit/msdyx/test/'
path = '/home/zezhou/McGillResearch/2019Manuscript_Analysis/femsimulation/iterationdata_sq/BMsimulation/Exponential/MSDfit/CalibratedFitting/' # Calibrated
os.chdir(path)
dataset = ['ecc0', 'ecc06', 'ecc08', 'ecc09', 'ecc095']

class msd_alpha:
    def __init__(self,name):
        self.name = name
        return

    def print(self):
        for item in dataset:
            file = np.loadtxt(self.name+item+'.txt')
            print('Alpha from clips: '+self.name+item)
            print(file)
        return

    def mean(self):
        for item in dataset:
            file = np.loadtxt(self.name+item+'.txt')
            mv = np.mean(file, axis=0)
            print('Mean of Alpha of '+ self.name+item+': '+str(mv))

    def std(self):
        for item in dataset:
            file = np.loadtxt(self.name+item+'.txt')
            mv = np.std(file, axis=0)
            print('Std of Alpha of '+ self.name+item+': '+str(mv))

    def sem(self):
        for item in dataset:
            file = np.loadtxt(self.name+item+'.txt')
            mv = np.std(file, axis=0)/np.sqrt(file.shape[0])
            print('Sem of Alpha of '+ self.name+item+': '+str(mv))

msdname = 'msdfit_clip_'
msdxname = 'msdfitx_clip_'
msdyname = 'msdfity_clip_'

a_msd = msd_alpha(name=msdname)
a_msd.sem()
a_msd.mean()

a_msdx = msd_alpha(name=msdxname)
a_msdx.sem()
a_msdx.mean()

a_msdy = msd_alpha(name=msdyname)
a_msdy.sem()
a_msdy.mean()