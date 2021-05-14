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

            
for j in range(len(first_row)):
    for file in os.listdir('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec'):
        #if file.startswith('1000000') and file.endswith('.dat'):
        #if file.startswith('1000010') and file.endswith('.dat'):
        if file.startswith('1001370') and file.endswith('.dat'):
            print(file)
            data=np.genfromtxt(os.path.join('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec', file), skip_footer=1)
            #values = open('E_1000000_%i.txt' %(j+1), 'w')
            #values = open('E_1000010_%i.txt' %(j+1), 'w')
            values = open('E_1001370_%i.txt' %(j+1), 'w')
            
            
            values.write(' #Energy Kinetic_Energy PotentiaL_Energy R \n')
            array_r=[]

            x=data[int(first_row[j]):int(last_row[j]),0] * ru
            y=data[int(first_row[j]):int(last_row[j]),1] * ru
            z=data[int(first_row[j]):int(last_row[j]),2] * ru
            vx=data[int(first_row[j]):int(last_row[j]),3] * vu
            vy=data[int(first_row[j]):int(last_row[j]),4] * vu
            vz=data[int(first_row[j]):int(last_row[j]),5] * vu
            m=data[int(first_row[j]):int(last_row[j]),6] * mu
            position=len(x)-1
                                        
            xcen=np.median(x)
            ycen=np.median(y)
            zcen=np.median(z)
            
            vxcen=np.median(vx)
            vycen=np.median(vy)
            vzcen=np.median(vz)
            
            print('Centre of mass',xcen,ycen,zcen)
            print('Sistemic velocity',vxcen,vycen,vzcen, np.abs(vxcen),np.abs(vycen),np.abs(vzcen))
            
            vx=vx-vxcen
            vy=vy-vycen
            vz=vz-vzcen
            #------------First particle-------------
            xpos=x[0]
            ypos=y[0]
            zpos=z[0]
            vxpos=vx[0]
            vypos=vy[0]
            vzpos=vz[0]
            mpos=m[0]
            
            #------------Other particle of the GC-------------
            x=x[1:]
            y=y[1:]
            z=z[1:]
            vx=vx[1:]
            vy=vy[1:]
            vz=vz[1:]
            m=m[1:]

            #Kinetic Energy
            K=cinetica(vxpos, vypos, vzpos)
            #Distances between first star and all the others
            vectorr= ((x-xpos)**2   +  (y-ypos)**2    +    (z-zpos)**2)**(0.5)
            array_r.append(vectorr)
            #Potential Energy
            U=potencial(G,m,vectorr)
            #Total Energy calculation
            E=K+U
            #Energy of the whole system
            ESYSTEM=K+(U/2)
            values.write(str(format(float(E) , '.4f'))  + ' ' + str(format(float(K) , '.4f'))  + ' ' + str(format(float(U) , '.10f'))  + ' ' + str(format(float(rcen(xpos,ypos,zpos,xcen,ycen,zcen)) , '.4f'))  + ' ' + '\n')
            
            for i in range(len(x)):
                xpos=x[0]
                ypos=y[0]
                zpos=z[0]
                vxpos=vx[0]
                vypos=vy[0]
                vzpos=vz[0]
                mpos=m[0]
                
                x=x[1:]
                y=y[1:]
                z=z[1:]
                vx=vx[1:]
                vy=vy[1:]
                vz=vz[1:]
                
                K=cinetica(vxpos, vypos, vzpos)
                vectorr=    ((x-xpos)**2   +     (y-ypos)**2    +         (z-zpos)**2)**(0.5)
                missing_r= [g[-position] for g in array_r]
                position=position-1
                array_r.append(vectorr)
                
                Upos=potencial(G, m[i+1:],vectorr)
                Uneg=potencial(G, m[:len(missing_r)],missing_r)
                
                E=K+Upos+Uneg
                ESYSTEM=ESYSTEM+K+(Upos+Uneg/2)
                values.write(str(format(float(E) , '.4f'))  + ' ' + str(format(float(K) , '.4f'))  + ' ' + str(format(float(Upos+Uneg) , '.10f'))  + ' ' + str(format(float(rcen(xpos,ypos,zpos,xcen,ycen,zcen)) , '.4f'))  + ' ' + '\n')
                
            print('Total Energy of the system')
            print(ESYSTEM)
