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
import os

def f(x, a, b, c):
    return a * np.exp(-(x - b)**2.0 / (2 * c**2))
    
first_row, last_row = np.loadtxt('/Volumes/Gabriela_Pink/GCproperties.dat', usecols=(5,6), unpack=True)
first_row=first_row+1

filenames= ['E_1001370_1.txt' , 'E_1001370_2.txt' , 'E_1001370_3.txt' , 'E_1001370_4.txt' , 'E_1001370_5.txt' , 'E_1001370_6.txt' , 'E_1001370_7.txt' , 'E_1001370_8.txt' , 'E_1001370_9.txt' , 'E_1001370_10.txt' , 'E_1001370_11.txt' , 'E_1001370_12.txt' , 'E_1001370_13.txt' , 'E_1001370_14.txt' , 'E_1001370_15.txt' , 'E_1001370_16.txt' , 'E_1001370_17.txt' , 'E_1001370_18.txt' , 'E_1001370_19.txt' , 'E_1001370_20.txt' , 'E_1001370_21.txt' , 'E_1001370_22.txt' , 'E_1001370_23.txt' , 'E_1001370_24.txt' , 'E_1001370_25.txt' , 'E_1001370_26.txt' , 'E_1001370_27.txt' , 'E_1001370_28.txt' , 'E_1001370_29.txt' , 'E_1001370_30.txt' , 'E_1001370_31.txt' , 'E_1001370_32.txt' , 'E_1001370_33.txt' , 'E_1001370_34.txt' , 'E_1001370_35.txt' , 'E_1001370_36.txt' , 'E_1001370_37.txt' , 'E_1001370_38.txt' , 'E_1001370_39.txt', 'E_1001370_40.txt', 'E_1001370_41.txt', 'E_1001370_42.txt']

aux_max=np.loadtxt('max_values_1.txt')
aux_max = aux_max.astype(int)
aux_mid=np.loadtxt('mid_values_5.txt')
aux_mid = aux_mid.astype(int)
aux_min=np.loadtxt('min_values_20.txt')
aux_min = aux_min.astype(int)
    
pdfEN1 = PdfPages('Energy_distribution_bins_check_s1370.pdf')

for i in range(len(first_row)):
    print(filenames[i])
    
    EN =np.loadtxt(filenames[i], usecols=(0), unpack=True)
    bins_ =np.linspace(min(EN), max(EN), 100)
    data_EN=plt.hist(EN, bins =bins_ )
    bines_EN=np.zeros(len(data_EN[0]))
    for j in range(len(bines_EN)):
        bines_EN[j]=(data_EN[1][j] + data_EN[1][j+1])/2
    
    aux1=np.where(aux_max>first_row[i])[0]
    aux2=np.where(aux_max<last_row[i])[0]
    valmax=list(set(aux1) & set(aux2))
    aux3=np.where(aux_mid>first_row[i])[0]
    aux4=np.where(aux_mid<last_row[i])[0]
    valmid=list(set(aux3) & set(aux4))
    aux5=np.where(aux_min>first_row[i])[0]
    aux6=np.where(aux_min<last_row[i])[0]
    valmin=list(set(aux5) & set(aux6))
    
    valmax=aux_max[valmax]-int(first_row[i])
    valmid=aux_mid[valmid]-int(first_row[i])
    valmin=aux_min[valmin]-int(first_row[i])

    fig=plt.figure(figsize=(10,10))
    plt.plot(bines_EN, data_EN[0], color = 'black' , marker='.')
    plt.hist(EN[valmax], bins =np.linspace(min(EN), max(EN), 100) , facecolor='white',histtype='stepfilled', edgecolor='red')
    plt.hist(EN[valmid], bins =np.linspace(min(EN), max(EN), 100) , facecolor='white',histtype='stepfilled', edgecolor='green')
    plt.hist(EN[valmin], bins =np.linspace(min(EN), max(EN), 100) , facecolor='white',histtype='stepfilled', edgecolor='blue')
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E', fontsize=16)
    txt="Caption: Total energy distribution of the cluster No %i in the snapshot number 1370. The blue, green and red histograms lines show show the current energy distribution of the subpopulations selected in the first snapshot for three set of energy ranges."%((i+1))
    plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig(pdfEN1, format='pdf')
    plt.close()

pdfEN1.close()
