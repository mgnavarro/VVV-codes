# -*- coding: utf-8 -*-/Users/gabrielanavarro/Dropbox
from pylab import *
import scipy
import numpy
import math
from math import sqrt
from decimal import *
from statistics import median
from scipy.optimize import curve_fit
from itertools import chain
from matplotlib.backends.backend_pdf import PdfPages
import random
from io import StringIO
import os

def rcen(x,y,z,x0,y0,z0):
    return sqrt(((x-x0)**2) + ((y-y0)**2) + ((z-z0)**2))

def cinetica(vx,vy,vz):
    return 0.5 * (vx**2+vy**2+vz**2)

def potencial(Gp,mass,dist):
    return -1*Gp * np.sum(np.divide(mass,dist))

#G=1
G=6.67 * (10**(-11))
mu=(10**6) * 1.99 * (10**(30))
ru=20 * 3.086 * (10**(16))
vu= (G * mu / ru)**0.5

first_row, last_row = np.loadtxt('/Volumes/Gabriela_Pink/GCproperties.dat', usecols=(5,6), unpack=True)
first_row=first_row+1

values = open('clusters_movements.txt', 'w')
#values.write(' #Energy Kinetic_Energy PotentiaL_Energy R \n')
            
            
for j in range(len(first_row)):
    for file in os.listdir('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec'):
        if file.startswith('100') and file.endswith('.dat'):
        #if file.startswith('100000') and file.endswith('.dat'):
            print('Snapshot ', file, ' for the cluster ', j, '...')
            data=np.genfromtxt(os.path.join('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec', file), skip_footer=1)

            array_r=[]

            x=data[int(first_row[j]):int(last_row[j]),0] * ru
            y=data[int(first_row[j]):int(last_row[j]),1] * ru
            z=data[int(first_row[j]):int(last_row[j]),2] * ru
            vx=data[int(first_row[j]):int(last_row[j]),3] * vu
            vy=data[int(first_row[j]):int(last_row[j]),4] * vu
            vz=data[int(first_row[j]):int(last_row[j]),5] * vu
            m=data[int(first_row[j]):int(last_row[j]),6] * mu
                                        
            xcen=np.median(x)
            ycen=np.median(y)
            zcen=np.median(z)
            
            vxcen=np.median(vx)
            vycen=np.median(vy)
            vzcen=np.median(vz)
            
            values.write(str(format(float(xcen)))  + ' ' + str(format(float(ycen)))  + ' ' + str(format(float(zcen))) + ' '  )
            values.write(str(format(float(vxcen) , '.4f'))  + ' ' + str(format(float(vycen) , '.4f'))  + ' ' + str(format(float(vzcen) , '.4f'))  + ' ')
            
    values.write( '\n')
