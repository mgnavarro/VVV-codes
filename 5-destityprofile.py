# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
from scipy.stats import sem
from scipy import stats
from scipy import interpolate
import statistics
import os

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
    
def create_circle(rad):
    circle= plt.Circle((0,0), radius= rad, fill=False, color='red')
    return circle

def show_shape(patch):
    ax=plt.gca()
    ax.add_patch(patch)
    plt.axis('scaled')
    plt.show()

            
first_row, last_row = np.loadtxt('/Volumes/Gabriela_Pink/GCproperties.dat', usecols=(5,6), unpack=True)
first_row=first_row+1

aux_max=np.loadtxt('max_values_20.txt')
aux_max = aux_max.astype(int)
aux_mid=np.loadtxt('mid_values_5.txt')
aux_mid = aux_mid.astype(int)
aux_min=np.loadtxt('min_values_1.txt')
aux_min = aux_min.astype(int)
#ds = open('Densityprof_vol_max_10.txt', 'w')
pdfEN1 = PdfPages('Volumetric_density_profile_bins.pdf')

G=6.67 * (10**(-11))
mu=(10**6) * 1.99 * (10**(30))
ru=20 * 3.086 * (10**(16))
vu= (G * mu / ru)**0.5

for file in os.listdir('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec'):
    #if file.startswith('1001369') and file.endswith('.dat'):
    if file.startswith('1000000') and file.endswith('.dat'):
        for i in range(len(first_row)):
        #for i in range(5):
            print('Computing the density profile for cluster number #: ', i, '...')
            data = np.genfromtxt(os.path.join('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec', file), skip_footer=1)
            
            '''
            x=data[int(first_row[j]):int(last_row[j]),0] * ru
            y=data[int(first_row[j]):int(last_row[j]),1] * ru
            z=data[int(first_row[j]):int(last_row[j]),2] * ru
            '''

            cmxp=np.mean(data[int(first_row[i]):int(last_row[i]),0])
            cmyp=np.mean(data[int(first_row[i]):int(last_row[i]),1])
            cmzp=np.mean(data[int(first_row[i]):int(last_row[i]),2])
            
            #COMPUTING THE CENTRE OF MASS OF THE CLUSTER IN PC
            cmx=(np.rad2deg(np.arctan((cmxp*20)/8178)))
            cmy=(np.rad2deg(np.arctan((cmyp*20)/8178)))
            cmz=(np.rad2deg(np.arctan((cmzp*20)/8178)))
            
            #Selecting the subpopulations of the cluster using the indexes
            aux1=np.where(aux_max>first_row[i])[0]
            aux2=np.where(aux_max<last_row[i])[0]
            valmax=list(set(aux1) & set(aux2))
            xmax=data[aux_max[valmax],0]
            ymax=data[aux_max[valmax],1]
            zmax=data[aux_max[valmax],2]
            
            aux3=np.where(aux_mid>first_row[i])[0]
            aux4=np.where(aux_mid<last_row[i])[0]
            valmid=list(set(aux3) & set(aux4))
            xmid=data[aux_mid[valmid],0]
            ymid=data[aux_mid[valmid],1]
            zmid=data[aux_mid[valmid],2]
            
            aux5=np.where(aux_min>first_row[i])[0]
            aux6=np.where(aux_min<last_row[i])[0]
            valmin=list(set(aux5) & set(aux6))
            xmin=data[aux_min[valmin],0]
            ymin=data[aux_min[valmin],1]
            zmin=data[aux_min[valmin],2]
            
            #All the positions in degrees for X Y and Z
            longmax=(np.rad2deg(np.arctan((xmax*20)/8178)))
            latmax=(np.rad2deg(np.arctan((ymax*20)/8178)))
            distancemax=(np.rad2deg(np.arctan((zmax*20)/8178)))
            rmax= (( (longmax-cmx)**2) + ((latmax-cmy)**2)  + ((distancemax-cmz)**2)     )**(0.5)
            
            longmid=(np.rad2deg(np.arctan((xmid*20)/8178)))
            latmid=(np.rad2deg(np.arctan((ymid*20)/8178)))
            distancemid=(np.rad2deg(np.arctan((zmid*20)/8178)))
            rmid= (((longmid-cmx)**2) + ((latmid-cmy)**2)  + ((distancemid-cmz)**2)     )**(0.5)
            
            longmin=(np.rad2deg(np.arctan((xmin*20)/8178)))
            latmin=(np.rad2deg(np.arctan((ymin*20)/8178)))
            distancemin=(np.rad2deg(np.arctan((zmin*20)/8178)))
            rmin= (((longmin-cmx)**2) + ((latmin-cmy)**2)  + ((distancemin-cmz)**2)     )**(0.5)
            
            
            areas=35
            r_plot=np.zeros(areas)
            RRLcorr_areamax=np.zeros(areas)
            RRLcorr_areamid=np.zeros(areas)
            RRLcorr_areamin=np.zeros(areas)
            rdef_Exp =np.linspace(np.log10((0.001)),np.log10((0.2)), areas)
            #rdef_Exp =np.linspace(np.log10(min(rmax)),np.log10(max(rmax)), areas)
            #ww_Exp=np.round(rdef_Exp)
            rdef=10**rdef_Exp
        
            for h in range(areas):
                rmin_=rdef[h]-(0.01/2)
                rmax_=rdef[h]+(0.01/2)
                r_plot[h]=rdef[h]
                vol_shell= ((4/3 * np.pi * (rmax_**3)) - (4/3 * np.pi * (rmin_**3)))

                aux1max=np.where(rmax>rmin_)[0]
                aux2max=np.where(rmax<rmax_)[0]
                valmax=list(set(aux1max) & set(aux2max))
                RRLcorr_areamax[h] = len(valmax) / vol_shell
                
                aux1mid=np.where(rmid>rmin_)[0]
                aux2mid=np.where(rmid<rmax_)[0]
                valmid=list(set(aux1mid) & set(aux2mid))
                RRLcorr_areamid[h] = len(valmid) / vol_shell
                
                aux1min=np.where(rmin>rmin_)[0]
                aux2min=np.where(rmin<rmax_)[0]
                valmin=list(set(aux1min) & set(aux2min))
                RRLcorr_areamin[h] = len(valmin) / vol_shell

            
            fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(22, 14))
            ax1.plot(longmax, latmax, color = 'red' ,  marker='o',linestyle=' ', markersize=1)
            ax1.plot(longmid, latmid, color = 'green' ,  marker='o',linestyle=' ', markersize=1)
            ax1.plot(longmin, latmin, color = 'blue' ,  marker='o',linestyle=' ', markersize=1)
            #plt.title(r'Cluster No %i -- #Part: %i' %((i+1),len(EN)), fontsize=16)
            ax1.set_ylabel(r'latitude (deg)', fontsize=16)
            ax1.set_xlabel(r'longitude (deg)', fontsize=16)
            ax1.set_box_aspect(1)
            ax1.set_xlim([np.mean(longmin)-0.22, np.mean(longmin)+0.22])
            ax1.set_ylim([np.mean(latmax)+0.22, np.mean(latmax)-0.22])
            #ax1.axis('off')
            ax1.text(0.5,-0.15, "Caption: Spatial distribution of cluster number %i in the first snapshot of the simulation. \n The blue, green and red distributions represent the most negative energies, mid energies and energies at the peak of the \n distribution, respectively."%((i+1)), size=12, ha="center", transform=ax1.transAxes)
            ax2.plot(r_plot, RRLcorr_areamax, color = 'red' ,  marker='o',linestyle=' ', markersize=8)
            ax2.plot(r_plot, RRLcorr_areamid, color = 'green' ,  marker='o',linestyle=' ', markersize=8)
            ax2.plot(r_plot, RRLcorr_areamin, color = 'blue' , marker='o',linestyle=' ', markersize=8)
            #plt.title(r'Cluster No %i -- #Part: %i' %((i+1),len(EN)), fontsize=16)
            ax2.set_ylabel(r'$\rho$', fontsize=16)
            ax2.set_xlabel(r'R (deg)', fontsize=16)
            ax2.set_box_aspect(1)
            ax2.set_yscale("log")
            ax2.set_xscale("log")
            ax2.set_xlim([0.0008, 0.2])
            ax2.set_ylim([100, 1000000000])
            #ax2.axis('off')
            ax2.text(0.5,-0.18, "Caption: Density profile for the three subpopulations of cluster number %i in the first snapshot of \n the simulation. The blue, green and red points represent the most negative energies, mid energies and energies in the \n peak of the distribution, respectively. The convertion from parsecs to degrees was implemented assuming the Galactocentric distance \n R0 = 8.178 ± 13stat. ± 22sys. pc. (The GRAVITY Collaboration, A&A 625, L10 (2019))."%((i+1)),  size=12, ha="center", transform=ax2.transAxes)
            fig.tight_layout()
            fig.savefig(pdfEN1, format='pdf')
            
pdfEN1.close()

'''
radius=np.zeros(len(RRLcorr_area))
radius=np.around(r_plot, decimals=2)
#radius = " ".join(str(v) for v in radius)
#radius = radius + '\n'
#rs.write(radius)
'''
