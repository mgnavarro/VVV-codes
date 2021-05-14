import pylab as plt
import numpy as np
import scipy
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from pylab import *
import scipy
import numpy
import math
from math import sqrt
from decimal import *
from statistics import median
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import numpy.polynomial.polynomial as poly
import matplotlib.mlab as mlab

def f(x, a, b, c):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2))

first_row, last_row = np.loadtxt('/Volumes/Gabriela_Pink/GCproperties.dat', usecols=(5,6), unpack=True)
first_row=first_row+1

cmap = plt.cm.get_cmap("hsv", len(first_row)+1)
pdf = PdfPages('Energy_distributions.pdf')

ex=np.zeros(len(first_row))
ey=np.zeros(len(first_row))
ux=np.zeros(len(first_row))
uy=np.zeros(len(first_row))
kx=np.zeros(len(first_row))
ky=np.zeros(len(first_row))
rx=np.zeros(len(first_row))
ry=np.zeros(len(first_row))

for i in range(len(first_row)):
#for i in range(2):
    EN, EK, U, R =np.loadtxt('E_1000000_%i.txt'%(i+1), usecols=(0,1,2,3), unpack=True)
            
    #---------------------------------EN---------------------------------
    VARIABLE=EN
    MIN, MAX = min(VARIABLE),max(VARIABLE)
    
    '''
    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()
    
    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.gca().set_xscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()
    '''
    
    fig = plt.figure(figsize=(22, 14))


    plt.subplot(2, 2, 1)
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E (J)', fontsize=16)
    plt.ylim([0.2, 800])
    plt.xlim([-10603579410.7111, 6272889821.5257])
    ex[i]=min(dataFQ[1])
    ey[i]=min(dataFQ[0])

    #plt.text(700, 0.5, 'Caption: Energy distribution for the cluster No %i. '%(i+1) , horizontalalignment='center',verticalalignment='center', size=15)


    #---------------------------------U---------------------------------
    VARIABLE=U
    MIN, MAX = min(VARIABLE),max(VARIABLE)
    
    '''
    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'U', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'U', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.gca().set_xscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'U', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()
    '''
    
    plt.subplot(2, 2, 2)
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'U (J)', fontsize=16)
    plt.ylim([1.0, 423.0])
    plt.xlim([-10810086680.387815, -49771455.482097924])
    ux[i]=min(dataFQ[1])
    uy[i]=min(dataFQ[0])
    #---------------------------------K---------------------------------

    VARIABLE=EK
    MIN, MAX = min(VARIABLE),max(VARIABLE)

    '''
    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'K', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50), facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_xscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'K', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50), facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_xscale("log")
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'K', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()
    '''
    
    plt.subplot(2, 2, 3)
    dataFQ=plt.hist(VARIABLE, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50), facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'K (J)', fontsize=16)
    plt.ylim([0.2, 2890.0])
    plt.xlim([169287.6504999999, 14041685904.957582])
    plt.text(200, 0.015, 'Caption: Energy distribution for cluster No %i. The top left and right panels show the total and potential energy distribution of the cluster, respectively. \n The lower left and right panels show the kinetic energy and radial distribution. All the distributions were computed using the International System of Units \n so the energies are presented in Joules.'%(i+1) , size=16)
    
    kx[i]=min(dataFQ[1])
    ky[i]=min(dataFQ[0])

    #---------------------------------R---------------------------------

    VARIABLE=R
    MIN, MAX = min(VARIABLE),max(VARIABLE)

    '''
    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins =np.linspace(MIN, MAX, 100) , facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'R', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50), facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_xscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'R', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()

    fig=plt.figure()
    dataFQ=plt.hist(VARIABLE, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50), facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_xscale("log")
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'R', fontsize=16)
    plt.title(r'Cluster No %i'%(i+1), fontsize=16)
    plt.savefig(pdf, format='pdf')
    plt.close()
    '''
    
    plt.subplot(2, 2, 4)
    dataFQ=plt.hist(VARIABLE, bins = 10 ** np.linspace(np.log10(MIN), np.log10(MAX), 50), facecolor='white',histtype='stepfilled', edgecolor='black')
    plt.gca().set_yscale("log")
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'R (m)', fontsize=16)
    plt.ylim([0.2, 1094.0])
    plt.xlim([215610782437246.1, 9.69183513908815e+17])
    rx[i]=min(dataFQ[1])
    ry[i]=min(dataFQ[0])
    
    
    plt.savefig(pdf, format='pdf')
    plt.close()
    

pdf.close()

'''
print(min(ex), min(ey))
print(min(ux), min(uy))
print(min(kx), min(ky))
print(min(rx), min(ry))
'''
