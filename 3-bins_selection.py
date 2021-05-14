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

min_porcent=1
mid_porcent=5
max_porcent=20

filenames= ['E_1000000_1.txt' , 'E_1000000_2.txt' , 'E_1000000_3.txt' , 'E_1000000_4.txt' , 'E_1000000_5.txt' , 'E_1000000_6.txt' , 'E_1000000_7.txt' , 'E_1000000_8.txt' , 'E_1000000_9.txt' , 'E_1000000_10.txt' , 'E_1000000_11.txt' , 'E_1000000_12.txt' , 'E_1000000_13.txt' , 'E_1000000_14.txt' , 'E_1000000_15.txt' , 'E_1000000_16.txt' , 'E_1000000_17.txt' , 'E_1000000_18.txt' , 'E_1000000_19.txt' , 'E_1000000_20.txt' , 'E_1000000_21.txt' , 'E_1000000_22.txt' , 'E_1000000_23.txt' , 'E_1000000_24.txt' , 'E_1000000_25.txt' , 'E_1000000_26.txt' , 'E_1000000_27.txt' , 'E_1000000_28.txt' , 'E_1000000_29.txt' , 'E_1000000_30.txt' , 'E_1000000_31.txt' , 'E_1000000_32.txt' , 'E_1000000_33.txt' , 'E_1000000_34.txt' , 'E_1000000_35.txt' , 'E_1000000_36.txt' , 'E_1000000_37.txt' , 'E_1000000_38.txt' , 'E_1000000_39.txt', 'E_1000000_40.txt', 'E_1000000_41.txt', 'E_1000000_42.txt']

min_values = open('min_values_%i.txt'%min_porcent, 'w')
mid_values = open('mid_values_%i.txt'%mid_porcent, 'w')
max_values = open('max_values_%i.txt'%max_porcent, 'w')

pdfEN1 = PdfPages('Energy_distribution_bins.pdf')

for i in range(len(first_row)):
#for i in range(2):
    print(filenames[i])
    
    EN =np.loadtxt(filenames[i], usecols=(0), unpack=True)
    bins_ =np.linspace(min(EN), max(EN), 100)
    data_EN=plt.hist(EN, bins =bins_ )
    bines_EN=np.zeros(len(data_EN[0]))
    for j in range(len(bines_EN)):
        bines_EN[j]=(data_EN[1][j] + data_EN[1][j+1])/2
    
    aux_pos=np.where(data_EN[0]==max(data_EN[0]))[0]
    aux_fit=np.where(bines_EN<=bines_EN[aux_pos][0])[0]
    z = np.polyfit(bines_EN[aux_fit], data_EN[0][aux_fit], 3)
    p = np.poly1d(z)
    axis_x =np.linspace(min(bines_EN[aux_fit]),max(bines_EN[aux_fit]), 400)
    
    auxlen_min=int(len(EN)*min_porcent/100)
    auxlen_mid=int(len(EN)*mid_porcent/100)
    auxlen_max=int(len(EN)*max_porcent/100)
    
    selected_aux=np.where(max(bines_EN[aux_fit])>=EN)[0]
    EN_ORDER=EN[selected_aux]
    EN_ORDER.sort()
    
    auxa=np.where(EN_ORDER[auxlen_min]>=EN)[0]
    aux3=np.where(EN_ORDER[  int( ( (len(EN_ORDER))/2) - (auxlen_mid/2) ) ]<=EN)[0]
    aux4=np.where(EN_ORDER[  int( ( (len(EN_ORDER))/2) + (auxlen_mid/2) )  ]>=EN)[0]
    auxb=list(set(aux3) & set(aux4))
    aux5=np.where(EN_ORDER[-(auxlen_max+1)]<=EN)[0]
    aux6=np.where(max(bines_EN[aux_fit])>=EN)[0]
    auxc=list(set(aux5) & set(aux6))
    
    '''
    fig=plt.figure()
    plt.plot(bines_EN, data_EN[0], color = 'black' , marker='.')
    plt.plot(bines_EN[aux_fit], data_EN[0][aux_fit], marker='.')
    plt.plot(axis_x, p(axis_x), color = 'red')
    plt.title(r'Cluster No %i -- #Part: %i' %((i+1),len(EN)), fontsize=16)
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E', fontsize=16)
    plt.savefig(pdfEN1, format='pdf')
    plt.close()
    '''
    
    fig=plt.figure(figsize=(10,10))
    plt.plot(bines_EN, data_EN[0], color = 'black' , marker='.')
    plt.axvline(x=EN_ORDER[0], color='blue')
    plt.axvline(x=EN_ORDER[auxlen_min], color='blue')
    plt.axvline(x=EN_ORDER[  int( ( (len(EN_ORDER))/2) - (auxlen_mid/2) ) ], color='green')
    plt.axvline(x=EN_ORDER[  int( ( (len(EN_ORDER))/2) + (auxlen_mid/2) )  ], color='green')
    plt.axvline(x=EN_ORDER[-auxlen_max], color='red')
    plt.axvline(x=EN_ORDER[-1], color='red')
    plt.hist(EN[auxa], bins =np.linspace(min(EN), max(EN), 100) , facecolor='white',histtype='stepfilled', edgecolor='blue')
    plt.hist(EN[auxb], bins =np.linspace(min(EN), max(EN), 100) , facecolor='white',histtype='stepfilled', edgecolor='green')
    plt.hist(EN[auxc], bins =np.linspace(min(EN), max(EN), 100) , facecolor='white',histtype='stepfilled', edgecolor='red')
    plt.ylabel(r'N', fontsize=16)
    plt.xlabel(r'E (J)', fontsize=16)
    plt.xlim([-10000000000, 5000000000])
    plt.ylim([1, 1000])

    plt.gca().set_yscale("log")
    txt="Caption: Total energy distribution for cluster No %i that contains %i particles. The blue, green and red vertical lines show the energy intervals for three sets of energy ranges corresponding to %i, %i and %i percent of the population in each case, respectively"%((i+1),len(EN), min_porcent, mid_porcent, max_porcent)
    plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig(pdfEN1, format='pdf')
    plt.close()
    

    for a in range(len(auxa)):
        min_values.write(str(format(int(auxa[a]+first_row[i]))) + ' ' + '\n')
        
    for b in range(len(auxb)):
        mid_values.write(str(format(int(auxb[b]+first_row[i]))) + ' ' + '\n')
        
    for c in range(len(auxc)):
        max_values.write(str(format(int(auxc[c]+first_row[i]))) + ' ' + '\n')

pdfEN1.close()
