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
    
def deg2rad(x):
    return 8.178*(np.tan(x * np.pi / 180))

def rad2deg(x):
    return (180*np.arctan(x/8.178) )/ np.pi
    
def best_fit(x, m, c):
    return 10**(m * np.log10(x) + c)
            
first_row, last_row = np.loadtxt('/Volumes/Gabriela_Pink/GCproperties.dat', usecols=(5,6), unpack=True)
first_row=first_row+1

aux_max=np.loadtxt('max_values_20.txt')
aux_max = aux_max.astype(int)
aux_mid=np.loadtxt('mid_values_5.txt')
aux_mid = aux_mid.astype(int)
aux_min=np.loadtxt('min_values_1.txt')
aux_min = aux_min.astype(int)

pdfEN1 = PdfPages('Projected_Density_profile_all.pdf')

for file in os.listdir('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec'):
    if file.startswith('1001369') and file.endswith('.dat'):
    #if file.startswith('1000000') and file.endswith('.dat'):
            data = np.genfromtxt(os.path.join('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec', file), skip_footer=1)
 
            x=data[0:int(last_row[-1]),0]
            y=data[0:int(last_row[-1]),1]
            
            #COMPUTING THE CENTRE OF MASS OF THE CLUSTER IN PC
            cmx=(np.rad2deg(np.arctan((np.mean(x)*20)/8178)))
            cmy=(np.rad2deg(np.arctan((np.mean(y)*20)/8178)))
                  
            #Selecting the 3 subpopulations using the indexes
            xmax=data[aux_max,0]
            ymax=data[aux_max,1]
            
            xmid=data[aux_mid,0]
            ymid=data[aux_mid,1]
            
            xmin=data[aux_min,0]
            ymin=data[aux_min,1]
            
            #All the positions in degrees for X and Y
            long=(np.rad2deg(np.arctan((x*20)/8178)))
            lat=(np.rad2deg(np.arctan((y*20)/8178)))
            r= (( (long-cmx)**2) + ((lat-cmy)**2)  )**(0.5)
                        
            longmax=(np.rad2deg(np.arctan((xmax*20)/8178)))
            latmax=(np.rad2deg(np.arctan((ymax*20)/8178)))
            rmax= (( (longmax-cmx)**2) + ((latmax-cmy)**2)  )**(0.5)
            
            longmid=(np.rad2deg(np.arctan((xmid*20)/8178)))
            latmid=(np.rad2deg(np.arctan((ymid*20)/8178)))
            rmid= (((longmid-cmx)**2) + ((latmid-cmy)**2) )**(0.5)
            
            longmin=(np.rad2deg(np.arctan((xmin*20)/8178)))
            latmin=(np.rad2deg(np.arctan((ymin*20)/8178)))
            rmin= (((longmin-cmx)**2) + ((latmin-cmy)**2) )**(0.5)
            
            areas=35
            r_plot=np.zeros(areas)
            RRLcorr_area=np.zeros(areas)
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
                area_shell= ((np.pi * (rmax_**2)) - (np.pi * (rmin_**2)))
                
                aux0a=np.where(r>rmin_)[0]
                aux0b=np.where(r<rmax_)[0]
                val=list(set(aux0a) & set(aux0b))
                RRLcorr_area[h] = len(val) / area_shell
                
                aux1max=np.where(rmax>rmin_)[0]
                aux2max=np.where(rmax<rmax_)[0]
                valmax=list(set(aux1max) & set(aux2max))
                RRLcorr_areamax[h] = len(valmax) / area_shell
                
                aux1mid=np.where(rmid>rmin_)[0]
                aux2mid=np.where(rmid<rmax_)[0]
                valmid=list(set(aux1mid) & set(aux2mid))
                RRLcorr_areamid[h] = len(valmid) / area_shell
                
                aux1min=np.where(rmin>rmin_)[0]
                aux2min=np.where(rmin<rmax_)[0]
                valmin=list(set(aux1min) & set(aux2min))
                RRLcorr_areamin[h] = len(valmin) / area_shell

            '''
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.scatter3D(data[0,0], data[0,1], data[0,2], c='red', marker='o', s=5)
            xa=data[0:int(first_row[cluster_name-1])-1,0]*20
            ya=data[0:int(first_row[cluster_name-1])-1,1]*20
            za=data[0:int(first_row[cluster_name-1])-1,2]*20
            xb=data[int(last_row[cluster_name-1])+1:478107,0]*20
            yb=data[int(last_row[cluster_name-1])+1:478107,1]*20
            zb=data[int(last_row[cluster_name-1])+1:478107,2]*20
            ax.scatter3D(xa, ya, za, c='black', marker='o', s=0.001)
            ax.scatter3D(xb, yb, zb, c='black', marker='o', s=0.001)
            ax.scatter3D(x, y, z, c='black', marker='o', s=0.001)
            ax.scatter3D(x[outer], y[outer], z[outer], c='blue', marker='o', s=0.05)
            ax.scatter3D(x[rrl], y[rrl], z[rrl], c='green', marker='o', s=0.05)
            ax.scatter3D(x[inner], y[inner], z[inner], c='red', marker='o', s=0.05)
            ax.set_xlabel('$X$ (pc)', fontsize=10)
            ax.set_ylabel('$Y$ (pc)', fontsize=10)
            ax.set_zlabel('$Z$ (pc)', fontsize=10, rotation = 0)
            ax.set_xlim3d(-1*cut,cut)
            ax.set_ylim3d(-1*cut,cut)
            ax.set_zlim3d(-1*cut,cut)
            #ax.yaxis._axinfo['label']['space_factor'] = 3.0
            plt.savefig(file[:-4] + '_' + str(cluster_name) +'.png')
            plt.close()

            fig=plt.figure(figsize=(12,10))
            plt.plot(long, lat, color = 'black' ,  marker='*',linestyle=' ', markersize=5)
            plt.plot(longmax, latmax, color = 'red' ,  marker='*',linestyle=' ', markersize=5)
            plt.plot(longmid, latmid, color = 'green' ,  marker='*',linestyle=' ', markersize=5)
            plt.plot(longmin, latmin, color = 'blue' ,  marker='*',linestyle=' ', markersize=5)
            #plt.title(r'Cluster No %i -- #Part: %i' %((i+1),len(EN)), fontsize=16)
            plt.xlabel(r'longitude (deg)', fontsize=16)
            plt.ylabel(r'latitude (deg)', fontsize=16)
            txt="Caption: Spatial distribution of all the 41 globular clusters in the last snapshot of the simulation. The blue, green and red represent the most negative energies, mid energies and energies in the peak of the distribution, respectively."
            plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
            plt.xlim([5, -5])
            plt.ylim([-5, 5])
            plt.savefig(pdfEN1, format='pdf')
            plt.close()
            '''

            fig=plt.figure(figsize=(14,12))
            plt.plot(r_plot, RRLcorr_area, color = 'black' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamax, color = 'red' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamid, color = 'green' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamin, color = 'blue' , marker='*',linestyle=' ', markersize=8)
            #plt.title(r'Cluster No %i -- #Part: %i' %((i+1),len(EN)), fontsize=16)
            plt.xlabel(r'$R \:(^\circ)$', fontsize=14)
            plt.ylabel(r'$\Sigma_{RRab} \: $ (per $deg^2$)', fontsize=14)
            plt.gca().set_yscale("log")
            plt.gca().set_xscale("log")
            txt="Caption: Projected density profile of all the 41 globular clusters in the last snapshot of the simulation. The black distribution correspond to the density \n profile for all the particles that belong to Globular Clusters. The blue, green and red represent the most negative energies, mid energies and energies in the \n peak of the distribution,  respectively. The conversion from parsecs to degrees was implemented assuming the Galactocentric \n distance R0 = 8.178 ± 13stat. ± 22sys. pc.  (The GRAVITY Collaboration, A&A 625, L10 (2019))."
            plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
            plt.savefig(pdfEN1, format='pdf')
            plt.close()

            #----------------------------------
            #FITTING THE DISTRIBUTION
            #----------------------------------
            cut= 0.01
            cutb=0.1
            
            aux_min=np.where(r_plot<cut)[0]
            m1=np.where(r_plot>cut)[0]
            m2=np.where(r_plot<cutb)[0]
            aux_mid=list(set(m1) & set(m2))
            aux_max=np.where(r_plot>cutb)[0]
            radiusint=r_plot[aux_min]
            radiusmid=r_plot[aux_mid]
            radiusext=r_plot[aux_max]
            
            axis_rrl1=np.linspace(min(radiusint), max(radiusint))
            axis_rrl2=np.linspace(min(radiusext), 3)
            axis_rrl3=np.linspace(min(radiusmid), max(radiusmid))

            dataint= RRLcorr_area[aux_min]
            datamid= RRLcorr_area[aux_mid]
            dataext= RRLcorr_area[aux_max]
            auxvalid=np.where(dataext>0)[0]
            p1, res1, _, _, _=np.polyfit(np.log10(radiusint), np.log10(dataint), 1, full=True)
            p2, res2, _, _, _=np.polyfit(np.log10(radiusext[auxvalid]), np.log10(dataext[auxvalid]), 1, full=True)
            p3, res3, _, _, _=np.polyfit(np.log10(radiusmid), np.log10(datamid), 1, full=True)
            
            dataint_max= RRLcorr_areamax[aux_min]
            datamid_max= RRLcorr_areamax[aux_mid]
            dataext_max= RRLcorr_areamax[aux_max]
            auxvalid=np.where(dataext_max>0)[0]
            p1max, res1, _, _, _=np.polyfit(np.log10(radiusint), np.log10(dataint_max), 1, full=True)
            p2max, res2, _, _, _=np.polyfit(np.log10(radiusext[auxvalid]), np.log10(dataext_max[auxvalid]), 1, full=True)
            p3max, res3, _, _, _=np.polyfit(np.log10(radiusmid), np.log10(datamid_max), 1, full=True)
            
            dataint_mid= RRLcorr_areamid[aux_min]
            datamid_mid= RRLcorr_areamid[aux_mid]
            dataext_mid= RRLcorr_areamid[aux_max]
            auxvalid=np.where(dataext_mid>0)[0]
            p1mid, res1, _, _, _=np.polyfit(np.log10(radiusint), np.log10(dataint_mid), 1, full=True)
            p2mid, res2, _, _, _=np.polyfit(np.log10(radiusext[auxvalid]), np.log10(dataext_mid[auxvalid]), 1, full=True)
            p3mid, res3, _, _, _=np.polyfit(np.log10(radiusmid), np.log10(datamid_mid), 1, full=True)
            
            dataint_min= RRLcorr_areamin[aux_min]
            datamid_min= RRLcorr_areamin[aux_mid]
            dataext_min= RRLcorr_areamin[aux_max]
            auxvalid=np.where(dataext_min>0)[0]
            p1min, res1, _, _, _=np.polyfit(np.log10(radiusint), np.log10(dataint_min), 1, full=True)
            p2min, res2, _, _, _=np.polyfit(np.log10(radiusext[auxvalid]), np.log10(dataext_min[auxvalid]), 1, full=True)
            p3min, res3, _, _, _=np.polyfit(np.log10(radiusmid), np.log10(datamid_min), 1, full=True)
        
            fig, ax1 = plt.subplots(constrained_layout=True, figsize=(14,12))
            ax1.loglog(r_plot, RRLcorr_area, marker='*', linestyle=' ',markersize=8, color='black')
            ax1.loglog(axis_rrl1, best_fit(axis_rrl1,p1[0],p1[1]), '-', color='black',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl2, best_fit(axis_rrl2,p2[0],p2[1]), '-', color='black',  lw=0.5,zorder=3)
            
            ax1.loglog(r_plot, RRLcorr_areamax, marker='*', linestyle=' ',markersize=8, color='red')
            ax1.loglog(axis_rrl1, best_fit(axis_rrl1,p1max[0],p1max[1]), '-', color='red',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl2, best_fit(axis_rrl2,p2max[0],p2max[1]), '-', color='red',  lw=0.5,zorder=3)
                        
            ax1.loglog(r_plot, RRLcorr_areamid, marker='*', linestyle=' ',markersize=8, color='green')
            ax1.loglog(axis_rrl1, best_fit(axis_rrl1,p1mid[0],p1mid[1]), '-', color='green',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl2, best_fit(axis_rrl2,p2mid[0],p2mid[1]), '-', color='green',  lw=0.5,zorder=3)
            
            ax1.loglog(r_plot, RRLcorr_areamin, marker='*', linestyle=' ',markersize=8, color='blue')
            ax1.loglog(axis_rrl1, best_fit(axis_rrl1,p1min[0],p1min[1]), '-', color='blue',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl2, best_fit(axis_rrl2,p2min[0],p2min[1]), '-', color='blue',  lw=0.5,zorder=3)
            
            ax1.loglog(axis_rrl3, best_fit(axis_rrl3,p3[0],p3[1]), '-', color='black',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl3, best_fit(axis_rrl3,p3max[0],p3max[1]), '-', color='red',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl3, best_fit(axis_rrl3,p3mid[0],p3mid[1]), '-', color='green',  lw=0.5,zorder=3)
            ax1.loglog(axis_rrl3, best_fit(axis_rrl3,p3min[0],p3min[1]), '-', color='blue',  lw=0.5,zorder=3)


            ax1.set_xlabel(r'$R \:(^\circ)$', fontsize=14)
            ax1.set_ylabel(r'$\Sigma_{RRab} \: $ (per $deg^2$)', fontsize=14)
            plt.text(0.003, 10**(2.7), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1[0]) , '.2f'), format((p3[0]) , '.2f'), format((p2[0]) , '.2f')), fontsize=17, color='black')
            
            plt.text(0.003, 10**(2.4), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1max[0]) , '.2f'), format((p3max[0]) , '.2f'), format((p2max[0]) , '.2f')), fontsize=17, color='red')
            
            plt.text(0.003, 10**(2.1), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1mid[0]) , '.2f'), format((p3mid[0]) , '.2f'), format((p2mid[0]) , '.2f')), fontsize=17, color='green')
                        
            plt.text(0.003, 10**(1.8), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1min[0]) , '.2f'), format((p3min[0]) , '.2f'), format((p2min[0]) , '.2f')), fontsize=17, color='blue')

            txt="Caption: Density profile of all the 41 globular clusters in the last snapshot of the simulation. The black distribution correspond to the \n density profile for all the  particles that belong to Globular Clusters. The blue, green and red represent the \n most negative energies, mid energies and energies in the peak of the distribution,  respectively. \n The lines correspond to the linear fit to the double power law in each case. The convertion from parsecs to degrees was implemented assuming \n the Galactocentric distance R0 = 8.178 ± 13stat. ± 22sys. pc.  (The GRAVITY Collaboration, A&A 625, L10 (2019))."
            plt.text(0.09, 0.1, txt, wrap=True, horizontalalignment='center', fontsize=12)
            
            secax = ax1.secondary_xaxis('top', functions=(deg2rad, rad2deg))
            secax.set_xlabel(r'$R \:$(kpc)', fontsize=14)
            plt.savefig(pdfEN1, format='pdf')
            plt.close()
            

            
            '''
            radius=np.zeros(len(RRLcorr_area))
            counts=np.zeros(len(RRLcorr_area))
            radius=np.around(r_plot, decimals=2)
            counts=np.around(RRLcorr_area, decimals=2)

            #radius = " ".join(str(v) for v in radius)
            #radius = radius + '\n'
            #rs.write(radius)
            
            counts = " ".join(str(v) for v in counts)
            counts = counts + '\n'
            ds.write(counts)
            '''

pdfEN1.close()
