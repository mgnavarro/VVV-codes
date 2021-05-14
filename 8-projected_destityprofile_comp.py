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

import matplotlib
from scipy import interpolate
from scipy import integrate
import math
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy import optimize
import decimal
from scipy.interpolate import interp1d
from matplotlib.transforms import Transform
from matplotlib.ticker import (AutoLocator, AutoMinorLocator)
from scipy.signal import find_peaks
from pylab import *
from scipy.optimize import curve_fit


def func_powerlaw(x, m, c, c0):
    return c0 + x**m * c

def func_powerlaw_area(x, m, c):
    return 10**(m * np.log10(x) + c)* (2 *np.pi * x)

def line_area(x, m, c):
    return (x*m + c ) * (2 *np.pi * x)

def solve_for_y(poly_coeffs, y):
    pc = poly_coeffs.copy()
    pc[-1] -= y
    return np.roots(pc)

def deg2rad(x):
    return 8.178*(np.tan(x * np.pi / 180))

def rad2deg(x):
    return (180*np.arctan(x/8.178) )/ np.pi

def tick_function(r_deg):
    r_rad= np.deg2rad(r_deg)
    aux_deg=np.tan(r_rad)
    rfinal=8.178* aux_deg
    return ["%.1f" % z for z in rfinal]

def tokpc(r_deg):
    r_rad= np.deg2rad(r_deg)
    aux_deg=np.tan(r_rad)
    rfinal=8.178* aux_deg
    return rfinal


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

r_tofit_rrl , counts_tofit_rrl ,err_rrl, err_rrl2, flag_rrl = np.loadtxt('table_deg.txt', usecols=(0,1,2,3,4), unpack=True)
colors=['white','red','blue','green','grey','purple','black']

sagitario_cut=23

first_cut=2.2
second_cut=8
third_cut=30

first_cutRC=2.2
second_cutRC=6.5
third_cutRC=30

RRAB35=5

r_tofit_rrl_cumulative=r_tofit_rrl
r_tofit_rrl=np.delete(r_tofit_rrl,0)
counts_tofit_rrl=np.delete(counts_tofit_rrl,0)
flag_rrl=np.delete(flag_rrl,0)
err_rrl=np.delete(err_rrl,0)

sgr1=np.where(r_tofit_rrl>sagitario_cut)[0]
sgr2=np.where(r_tofit_rrl<third_cut)[0]
sgr=list(set(sgr1) & set(sgr2))
r_tofit_rrl_sgr=r_tofit_rrl[sgr]
counts_tofit_rrl_sgr=counts_tofit_rrl[sgr]
flag_rrl_sgr=flag_rrl[sgr]
err_rrl_sgr=err_rrl[sgr]
nosgr1=np.where(r_tofit_rrl<sagitario_cut)[0]
nosgr2=np.where(r_tofit_rrl>third_cut)[0]
nosgr=list(set(nosgr1) | set(nosgr2))
r_tofit_rrl_nosgr=r_tofit_rrl[nosgr]
counts_tofit_rrl_nosgr=counts_tofit_rrl[nosgr]
flag_rrl_nosgr=flag_rrl[nosgr]
err_rrl_nosgr=err_rrl[nosgr]

print('---------------------------------------RRL---------------------------------------')
r1=np.where(r_tofit_rrl<first_cut)[0]
r2=np.where(r_tofit_rrl>first_cut)[0]
r3=np.where(r_tofit_rrl<second_cut)[0]
r4=list(set(r2) & set(r3))
r5=np.where(r_tofit_rrl>second_cut)[0]
r6=np.where(r_tofit_rrl<sagitario_cut)[0]
r7=list(set(r5) & set(r6))

axis_rrl1=np.linspace(1, first_cut+0.3)
axis_rrl2=np.linspace(first_cut-0.3, second_cut+1)
axis_rrl3=np.linspace(second_cut-1,third_cut)

axis_rrl1=np.linspace(1, first_cut)
axis_rrl2=np.linspace(first_cut, second_cut)
axis_rrl3=np.linspace(second_cut,third_cut)

#--------------------------FITS
p, res, _, _, _=np.polyfit(r_tofit_rrl[r1], counts_tofit_rrl[r1], 1, full=True)
p_plot = np.poly1d(p)
print('First Section linear fit',(p_plot))
print('First Section ',np.poly1d(np.polyfit(np.log10(r_tofit_rrl[r1]), np.log10(counts_tofit_rrl[r1]), 1)))
print('Second Section ',np.poly1d(np.polyfit(np.log10(r_tofit_rrl[r4]), np.log10(counts_tofit_rrl[r4]), 1)))
print('Third Section ',np.poly1d(np.polyfit(np.log10(r_tofit_rrl[r7]), np.log10(counts_tofit_rrl[r7]), 1)))

#--------------------------ERRORES
p1, res1, _, _, _=np.polyfit(np.log10(r_tofit_rrl[r1]), np.log10(counts_tofit_rrl[r1]), 1, full=True)
p2, res2, _, _, _=np.polyfit(np.log10(r_tofit_rrl[r4]), np.log10(counts_tofit_rrl[r4]), 1, full=True)
p3, res3, _, _, _=np.polyfit(np.log10(r_tofit_rrl[r7]), np.log10(counts_tofit_rrl[r7]), 1, full=True)

p1_=np.polyfit(np.log10(r_tofit_rrl[r1]), np.log10(counts_tofit_rrl[r1]), 1, cov=True)
p2_=np.polyfit(np.log10(r_tofit_rrl[r4]), np.log10(counts_tofit_rrl[r4]), 1, cov=True)
p3_=np.polyfit(np.log10(r_tofit_rrl[r7]), np.log10(counts_tofit_rrl[r7]), 1, cov=True)

popt1, pcov1 = curve_fit(func_powerlaw, r_tofit_rrl[r1], counts_tofit_rrl[r1],p0 = np.asarray([p1[0],p1[1],1000]))
error1 = []
for i in range(len(popt1)):
    try:
        error1.append(np.absolute(pcov1[i][i])**0.5)
    except:
        error1.append( 0.00 )
perr_curvefit1 = np.array(error1)

popt2, pcov2 = curve_fit(func_powerlaw, r_tofit_rrl[r2], counts_tofit_rrl[r2],p0 = np.asarray([p2[0],p2[1],1000]))
error2 = []
for i in range(len(popt2)):
    try:
        error2.append(np.absolute(pcov2[i][i])**0.5)
    except:
        error2.append( 0.00 )
perr_curvefit2 = np.array(error2)

popt3, pcov3 = curve_fit(func_powerlaw, r_tofit_rrl[r3], counts_tofit_rrl[r3],p0 = np.asarray([p3[0],p3[1],1000]))
error3 = []
for i in range(len(popt3)):
    try:
        error3.append(np.absolute(pcov3[i][i])**0.5)
    except:
        error3.append( 0.00 )
perr_curvefit3 = np.array(error3)

#print('E_First Section linear fit',res)
#print('E_First Section ',res1, sqrt(p1_[1][0][0]), perr_curvefit1)
#print('E_Second Section ',res2, sqrt(p2_[1][0][0]), perr_curvefit2)
#print('E_Third Section',res3, sqrt(p3_[1][0][0]), perr_curvefit3)
#print('\n')

#--------------------------COMPUTATION
#USING THE STRAIGTH LINE TO COMPUTE THE MIDDLE VALUE FROM THE MAX.
yaxis_mid_normal=np.poly1d(np.polyfit(r_tofit_rrl[r1], counts_tofit_rrl[r1], 1))(np.unique(0)) /2
r0_normal=solve_for_y(np.polyfit(r_tofit_rrl[r1], counts_tofit_rrl[r1], 1), yaxis_mid_normal)
print('Max nstars/area SL= ',np.poly1d(np.polyfit(r_tofit_rrl[r1], counts_tofit_rrl[r1], 1))(np.unique(0)) )
print('FHWM RRL SL(deg)', r0_normal)
print('FHWM RRL SL(kpc)', tokpc(r0_normal))

aaa, bbbb = quad(func_powerlaw_area, 0, 0.7, args=(p1[0],p1[1]))

#USING THE POWER LAW LINE TO COMPUTE THE MIDDLE VALUE FROM THE MAX.
yaxis_mid_normal_PL=best_fit(0.1,p1[0],p1[1]) /2
#print(best_fit(0.313,p1[0],p1[1]))
r0_normal_PL=0.313
print('Max nstars/area PL= ',best_fit(0.1,p1[0],p1[1]))
print('FHWM RRL PL(deg)', r0_normal_PL)
print('FHWM RRL PL(kpc)', tokpc(r0_normal_PL))

#INTEGRATING NUMBER OF STARS
ansl, errl = quad(line_area, 0, first_cut, args=(p[0],p[1]))
ans1, err1 = quad(func_powerlaw_area, 0, first_cut, args=(p1[0],p1[1]))
ans2, err2 = quad(func_powerlaw_area, first_cut, second_cut, args=(p2[0],p2[1]))
ans3, err3 = quad(func_powerlaw_area, second_cut, third_cut, args=(p3[0],p3[1]))
print('Total integral (SL)= ', ansl+ans2+ans3)
print('Total integral (PL)= ', ans1+ans2+ans3 )
print('\n')

ansnsceff, errnsceff = quad(line_area, 0, 0.03, args=(p[0],p[1]))
ansnsceffa, errnsceffa = quad(line_area, 0, 0.05, args=(p[0],p[1]))
ansnsc, errnsc = quad(line_area, 0, 0.25, args=(p[0],p[1]))
ansnb, errnb = quad(line_area, 0, 0.7, args=(p[0],p[1]))
ansn_rsum, errn_rsum = quad(line_area, 0, 1.7592, args=(p[0],p[1]))


print('------------------------------------------')
print('------------------------------------------')
print('------------------------------------------ \n \n \n')

print('Straight line')
print('Expected RRL within the effective radius of the NSC (4 PC) = ', ansnsceff ,'+-', errnsceff )
print('Expected RRL within the effective radius of the NSC (7 PC) = ', ansnsceffa ,'+-', errnsceffa )
print('Expected RRL in the NSC (35 PC) = ', ansnsc ,'+-', errnsc )
print('Porcentage (DONG E MINNITI)= ', RRAB35*100/ansnsc )
print('Expected RRL in the Nuclear Bulge (100 PC) = ',ansnb , '+-', errnb)
print('Density of RRL at galactocentric distance of 1.6= ',(np.poly1d(p_plot)(1.6)))
print('Expected RRL inside R (1.7592) = ', ansn_rsum ,'+-', errn_rsum )

ansnsceff, errnsceff = quad(func_powerlaw_area, 0, 0.03, args=(p1[0],p1[1]))
ansnsceffa, errnsceffa = quad(func_powerlaw_area, 0, 0.05, args=(p1[0],p1[1]))
ansnsc, errnsc = quad(func_powerlaw_area, 0, 0.25, args=(p1[0],p1[1]))
ansnb, errnb = quad(func_powerlaw_area, 0, 0.7, args=(p1[0],p1[1]))
ansn_rsum, errn_rsum = quad(func_powerlaw_area, 0, 1.7592, args=(p1[0],p1[1]))

print('Power law')
print('Expected RRL within the effective radius of the NSC (4 PC) = ', ansnsceff ,'+-', errnsceff )
print('Expected RRL within the effective radius of the NSC (7 PC) = ', ansnsceffa ,'+-', errnsceffa )
print('Expected RRL in the NSC (35 PC) = ', ansnsc ,'+-', errnsc )
print('Porcentage (DONG E MINNITI)= ', RRAB35*100/ansnsc )
print('Expected RRL in the Nuclear Bulge (100 PC) = ',ansnb , '+-', errnb)
print('Density of RRL at galactocentric distance of 1.6 = ',(best_fit(1.6,p1[0],p1[1])))
print('Expected RRL inside R (1.7592) = ', ansn_rsum ,'+-', errn_rsum )


print('------------------------------------------')
print('------------------------------------------')
print('------------------------------------------ \n \n \n')

forcomp, forcomperr = quad(func_powerlaw_area, 0, 2, args=(p1[0],p1[1]))
print('Expected RRL within 2 deg = ', forcomp ,'+-', forcomperr )


#--------------------------CUMULATIVE DISTRIBUTION
cum_rrl=np.zeros(len(r_tofit_rrl_cumulative))
cum_rrle=np.zeros(len(r_tofit_rrl_cumulative))
cum_rrl_counts=np.zeros(len(r_tofit_rrl_cumulative))
cum_rrl_counts_err=np.zeros(len(r_tofit_rrl_cumulative))
cum_flag_rrl=np.zeros(len(r_tofit_rrl_cumulative))

for i in range(0,len(r_tofit_rrl_cumulative)):
    if i == 0:
        cum_rrl[i]=0
        cum_rrle[i]=0
        cum_rrl_counts[i]= 0
        cum_rrl_counts_err[i]= 0
        cum_flag_rrl[i]=1
    else:
        if r_tofit_rrl_cumulative[i]<=first_cut:
            cum_rrl[i], cum_rrle[i] = quad(func_powerlaw_area, 0, r_tofit_rrl_cumulative[i], args=(p1[0],p1[1]))
            c_new= best_fit(r_tofit_rrl_cumulative[i],p1[0],p1[1])
        elif r_tofit_rrl_cumulative[i]>first_cut and r_tofit_rrl_cumulative[i]<=second_cut:
            prev_integral, prev_err=quad(func_powerlaw_area, 0, first_cut, args=(p1[0],p1[1]))
            new_value, new_err=quad(func_powerlaw_area, first_cut, r_tofit_rrl_cumulative[i], args=(p2[0],p2[1]))
            cum_rrl[i]=new_value + prev_integral
            cum_rrle[i]= prev_err+new_err
            c_new=best_fit(r_tofit_rrl_cumulative[i],p2[0],p2[1])
        else:
            prev_integral, prev_err=quad(func_powerlaw_area, 0, first_cut, args=(p1[0],p1[1]))
            sec_integral, sec_err =quad(func_powerlaw_area, first_cut, second_cut, args=(p2[0],p2[1]))
            new_value, new_err=quad(func_powerlaw_area, second_cut,r_tofit_rrl_cumulative[i], args=(p3[0],p3[1]))
            cum_rrl[i]=new_value + prev_integral + sec_integral
            cum_rrle[i]=  prev_err + new_err + sec_err
            c_new=best_fit(r_tofit_rrl_cumulative[i],p3[0],p3[1])
        cum_flag_rrl[i]=flag_rrl[i-1]
        cum_rrl_counts[i]= (((np.pi *r_tofit_rrl_cumulative[i]* r_tofit_rrl_cumulative[i]) - (np.pi * r_tofit_rrl_cumulative[i-1]* r_tofit_rrl_cumulative[i-1]) ) * c_new) + cum_rrl_counts[i-1]
        cum_rrl_counts_err[i]= (   ( ( ( (np.pi *r_tofit_rrl_cumulative[i]*r_tofit_rrl_cumulative[i]) - (np.pi*r_tofit_rrl_cumulative[i-1]* r_tofit_rrl_cumulative[i-1]) ) **2 ) * err_rrl[i-1]**2 ) + cum_rrl_counts_err[i-1]**2    ) **0.5

pdfEN1 = PdfPages('Projected_Density_profile_comparison.pdf')

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#------------------------------------------------PLOTS-----------------------------------------------
#-----------------------------------------------DENSITY----------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

#------------------------------------DENSITY PROFILE RRL
new_tick_locations = np.array([0, 5, 10, 15, 20,25, 30])
new_tick_locationsy = np.array([0, 100, 200, 300, 400, 500, 600,700 ])

new_tick_locations_log = np.array([10, 100])
fig = plt.figure(constrained_layout=True,figsize=(14,12))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
for i in range(len(r_tofit_rrl_nosgr)):
    ax1.errorbar(r_tofit_rrl_nosgr[i], counts_tofit_rrl_nosgr[i], yerr=err_rrl_nosgr[i], marker='.',markersize=17, color=colors[int(flag_rrl_nosgr[i])],capsize=1.6, elinewidth=0.6, alpha=0.4,zorder=1)
for i in range(len(r_tofit_rrl_sgr)):
    ax1.errorbar(r_tofit_rrl_sgr[i], counts_tofit_rrl_sgr[i], yerr=err_rrl_sgr[i], marker='o', mfc='none',markersize=10, color=colors[int(flag_rrl_sgr[i])],capsize=1.6, elinewidth=0.6 ,zorder=2)

ax1.plot(axis_rrl1, best_fit(axis_rrl1,p1[0],p1[1]), '-', color='black', lw=2,zorder=3)
ax1.plot(axis_rrl2, best_fit(axis_rrl2,p2[0],p2[1]),'g--', color='black', lw=2,zorder=4)
ax1.plot(axis_rrl3, best_fit(axis_rrl3,p3[0],p3[1]),  ':', color='black', lw=2,zorder=5)
ax1.set_xlabel(r'$R \:(^\circ)$', fontsize=14)
ax1.set_ylabel(r'$\Sigma_{RRab} \: $ (per sq.deg.)', fontsize=14)
ax1.set_ylim(0,max(counts_tofit_rrl)+ 50 )
ax1.set_xlim(0,max(r_tofit_rrl)+1)
ax1.set_xticks(new_tick_locations)
ax1.set_yticks(new_tick_locationsy)
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$R \:$(kpc)', fontsize=14)
plt.savefig(pdfEN1, format='pdf')
plt.close()


print('\n \n')
print('---------------------------------------FROM CUMULATIVE---------------------------------------')

print('---------------------------------------RRL---------------------------------------')

fig, ax1 = plt.subplots(constrained_layout=True, figsize=(6, 5))
for i in range(len(r_tofit_rrl_cumulative)):
    ax1.plot(r_tofit_rrl_cumulative[i],cum_rrl_counts[i], '.', color=colors[int(cum_flag_rrl[i])], markersize=8)
    ax1.plot(r_tofit_rrl_cumulative[i],cum_rrl[i], '*', color=colors[int(cum_flag_rrl[i])], markersize=8)
ax1.set_xticks(new_tick_locations)
plt.savefig(pdfEN1, format='pdf')
plt.close()

print('Number of stars at 30 deg:')
edge=np.where(r_tofit_rrl_cumulative>27)[0]
for30=np.polyfit(r_tofit_rrl_cumulative[edge], cum_rrl_counts[edge], 1)
p_for30 = np.poly1d(for30)
print('Linear fit:',p_for30(30))
totalstars=p_for30(30)
print('Polynomial fit:',np.poly1d(np.polyfit(r_tofit_rrl_cumulative, cum_rrl_counts, 3))(30))
print('Number of RRab for GC = ', p_for30(30)/86)
axis_rrl1=np.linspace(0, 35)

fig, ax1 = plt.subplots(constrained_layout=True, figsize=(6, 5))
ax1.plot(axis_rrl1,np.poly1d(np.polyfit(r_tofit_rrl_cumulative, cum_rrl_counts, 3))(axis_rrl1), markersize=6, color='blue')
ax1.plot(axis_rrl1,p_for30(axis_rrl1), color='red', markersize=6)
ax1.plot(r_tofit_rrl_cumulative, cum_rrl_counts, '*',markersize=6, color='black')
plt.savefig(pdfEN1, format='pdf')
plt.close()

print('Radius for half:')
f = interp1d(r_tofit_rrl_cumulative, cum_rrl_counts)
del1=np.where(cum_rrl_counts/totalstars<0.6)[0]
del2=np.where(cum_rrl_counts/totalstars>0.4)[0]
del3=list(set(del1) & set(del2))
r0cum=solve_for_y(np.polyfit(r_tofit_rrl_cumulative[del3], cum_rrl_counts[del3], 1)/totalstars, 0.5)
print(r0cum, tokpc(r0cum), f(r0cum)/totalstars)
print('Number of stars at 10 deg:')
print('Interpolation:', f(10))

for10=np.polyfit(r_tofit_rrl_cumulative[del3], cum_rrl_counts[del3], 1)
p_for10 = np.poly1d(for10)
print('Linear fit:',p_for10(10))

fig, ax1 = plt.subplots(constrained_layout=True, figsize=(6, 5))
ax1.plot(r_tofit_rrl_cumulative, cum_rrl_counts, '*',markersize=6, color='black')
ax1.plot(r_tofit_rrl_cumulative[del3], cum_rrl_counts[del3], '*',markersize=6, color='blue')
ax1.plot(axis_rrl1,p_for10(axis_rrl1), color='red', markersize=6)
plt.savefig(pdfEN1, format='pdf')
plt.close()


#------------------------------------CUMULATIVE PROFILE RRL
f = interp1d(r_tofit_rrl_cumulative, cum_rrl_counts)
fig, ax1 = plt.subplots(constrained_layout=True, figsize=(6, 5))
ax2 = ax1.twiny()
for i in range(len(r_tofit_rrl_cumulative)-1):
    ax1.errorbar(r_tofit_rrl_cumulative[i], cum_rrl_counts[i], yerr=cum_rrl_counts_err[i], marker='.', linestyle=' ', color=colors[int(cum_flag_rrl[i])],capsize=1.6, elinewidth=0.6 , markersize=6)
ax1.set_xlim(0,max(r_tofit_rrl)+1)
ax1.set_xlabel(r'$R \:(^\circ)$', fontsize=14)
ax1.set_ylabel(r'$N_{RRab}(<r)/N_{Total}$', fontsize=14)
ax1.set_xticks(new_tick_locations)
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$R \:$(kpc)', fontsize=14)
plt.savefig(pdfEN1, format='pdf')
plt.close()

fig, ax1 = plt.subplots(constrained_layout=True, figsize=(6, 5))
ax2 = ax1.twiny()
for i in range(len(r_tofit_rrl_cumulative)-1):
    ax1.errorbar(r_tofit_rrl_cumulative[i], cum_rrl_counts[i]/max(cum_rrl_counts), yerr=cum_rrl_counts_err[i]/max(cum_rrl_counts), marker='.', linestyle=' ', color=colors[int(cum_flag_rrl[i])],capsize=1.6, elinewidth=0.6 , markersize=14, alpha=0.7)
ax1.set_ylim(0,1.05)
ax1.set_xlim(0,max(r_tofit_rrl)+1)
ax1.set_xlabel(r'$R \:(^\circ)$', fontsize=14)
ax1.set_ylabel(r'$N_{RRab}(<r)/N_{Total}$', fontsize=14)
ax1.set_xticks(new_tick_locations)
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))
ax2.set_xlabel(r'$R \:$(kpc)', fontsize=14)
#plt.savefig('fig4.png')
plt.savefig(pdfEN1, format='pdf')
plt.close()



aux_max=np.loadtxt('max_values_20.txt')
aux_max = aux_max.astype(int)
aux_mid=np.loadtxt('mid_values_5.txt')
aux_mid = aux_mid.astype(int)
aux_min=np.loadtxt('min_values_1.txt')
aux_min = aux_min.astype(int)

for file in os.listdir('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec'):
    if file.startswith('1001369') and file.endswith('.dat'):
    #if file.startswith('1000000') and file.endswith('.dat'):
            data = np.genfromtxt(os.path.join('/Volumes/Gabriela_Pink/Data_FIRST_MEGANMW/exec', file), skip_footer=1)
 
            x=data[0:int(last_row[-1]),0]
            y=data[0:int(last_row[-1]),1]
            
            cmx=(np.rad2deg(np.arctan((np.mean(x)*20)/8178)))
            cmy=(np.rad2deg(np.arctan((np.mean(y)*20)/8178)))
                        
            xmax=data[aux_max,0]
            ymax=data[aux_max,1]
            
            xmid=data[aux_mid,0]
            ymid=data[aux_mid,1]
            
            xmin=data[aux_min,0]
            ymin=data[aux_min,1]
            
            
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
            
            areas=30
            r_plot=np.zeros(areas)
            RRLcorr_area=np.zeros(areas)
            RRLcorr_areamax=np.zeros(areas)
            RRLcorr_areamid=np.zeros(areas)
            RRLcorr_areamin=np.zeros(areas)
            rdef_Exp =np.linspace(np.log10(min(rmax)),np.log10(max(rmax)), areas)
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

        
            #------------------------------------DENSITY PROFILE RRL
            fig, ax1 = plt.subplots(constrained_layout=True, figsize=(14,12))
            for i in range(len(r_tofit_rrl_nosgr)):
                ax1.errorbar(r_tofit_rrl_nosgr[i], counts_tofit_rrl_nosgr[i], yerr=err_rrl_nosgr[i], marker='.',markersize=17, color=colors[int(flag_rrl_nosgr[i])],capsize=1.6, elinewidth=0.6, alpha=0.4,zorder=1)
            for i in range(len(r_tofit_rrl_sgr)):
                ax1.errorbar(r_tofit_rrl_sgr[i], counts_tofit_rrl_sgr[i], yerr=err_rrl_sgr[i], marker='o', mfc='none',markersize=10, color=colors[int(flag_rrl_sgr[i])],capsize=1.6, elinewidth=0.6, zorder=2 )
            '''
            for i in range(len(r_sim)):
                ax1.errorbar(r_sim[i], counts_sim[i], yerr=err_sim[i], marker='.', markersize=17, color=colors[int(flag_sim[i])],capsize=1.6, elinewidth=0.6, alpha=0.4,zorder=1)
            '''
            
            plt.plot(r_plot, RRLcorr_area, color = 'black' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamax, color = 'red' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamid, color = 'green' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamin, color = 'blue' , marker='*',linestyle=' ', markersize=8)
            
            ax1.loglog(axis_rrl1, best_fit(axis_rrl1,p1[0],p1[1]), '-', color='black',  lw=2,zorder=3)
            ax1.loglog(axis_rrl2, best_fit(axis_rrl2,p2[0],p2[1]),'g--', color='black', lw=2,zorder=4)
            ax1.loglog(axis_rrl3, best_fit(axis_rrl3,p3[0],p3[1]),  ':', color='black', lw=2,zorder=5)
            ax1.set_xlabel(r'$R \:(^\circ)$', fontsize=14)
            ax1.set_ylabel(r'$\Sigma_{RRab} \: $ (per sq.deg.)', fontsize=14)
            ax2.set_xlim(ax1.get_xlim())
            ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            secax = ax1.secondary_xaxis('top', functions=(deg2rad, rad2deg))
            secax.set_xlabel(r'$R \:$(kpc)', fontsize=14)
            txt="Caption: Projected density profile of all the 41 globular clusters in the last snapshot of the simulation. The black distribution correspond to the density \n profile for all the particles that belong to Globular Clusters. The blue, green and red represent the most negative energies, mid energies and energies in the \n peak of the distribution,  respectively. The conversion from parsecs to degrees was implemented assuming the Galactocentric \n distance R0 = 8.178 ± 13stat. ± 22sys. pc.  (The GRAVITY Collaboration, A&A 625, L10 (2019))."
            #plt.figtext(0.5, 0.000001, txt, wrap=True, horizontalalignment='center', fontsize=12)
            plt.savefig(pdfEN1, format='pdf')
            plt.close()

           
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
            axis_rrl2=np.linspace(min(radiusext), 7)
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
        
        
            
            ra, ea = quad(func_powerlaw_area, 0, 0.015, args=(p1min[0],p1min[1]))
            rb, eb = quad(func_powerlaw_area, 0.08, 2, args=(p2min[0],p2min[1]))
            rc, ec = quad(func_powerlaw_area, 0.015, 0.08, args=(p3min[0],p3min[1]))
            
            print('Expected RRL from simulation within 2 deg = ', ra+rb+rc ,'+-', ea+eb+ec )
            
            

            
            
            
            
            
            
            print('We need = ', (forcomp)/(ra+rb+rc)*41 , ' Globular clusters' )
            
            
            
            
            
            
            
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

            for i in range(len(r_tofit_rrl_nosgr)):
                ax1.errorbar(r_tofit_rrl_nosgr[i], counts_tofit_rrl_nosgr[i], yerr=err_rrl_nosgr[i], marker='.',markersize=17, color=colors[int(flag_rrl_nosgr[i])],capsize=1.6, elinewidth=0.6, alpha=0.4,zorder=1)
            for i in range(len(r_tofit_rrl_sgr)):
                ax1.errorbar(r_tofit_rrl_sgr[i], counts_tofit_rrl_sgr[i], yerr=err_rrl_sgr[i], marker='o', mfc='none',markersize=10, color=colors[int(flag_rrl_sgr[i])],capsize=1.6, elinewidth=0.6, zorder=2 )
            
            plt.plot(r_plot, RRLcorr_area, color = 'black' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamax, color = 'red' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamid, color = 'green' ,  marker='*',linestyle=' ', markersize=8)
            plt.plot(r_plot, RRLcorr_areamin, color = 'blue' , marker='*',linestyle=' ', markersize=8)
            
            #ax1.loglog(axis_rrl1, best_fit(axis_rrl1,p1[0],p1[1]), '-', color='black',  lw=2,zorder=3)
            #ax1.loglog(axis_rrl2, best_fit(axis_rrl2,p2[0],p2[1]),'g--', color='black', lw=2,zorder=4)
            #ax1.loglog(axis_rrl3, best_fit(axis_rrl3,p3[0],p3[1]),  ':', color='black', lw=2,zorder=5)
            
            ax1.set_xlabel(r'$R \:(^\circ)$', fontsize=14)
            ax1.set_ylabel(r'$\Sigma_{RRab} \: $ (per $deg^2$)', fontsize=14)
            plt.text(0.0003, 10**(2.7-3), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1[0]) , '.2f'), format((p3[0]) , '.2f'), format((p2[0]) , '.2f')), fontsize=17, color='black')
            
            plt.text(0.0003, 10**(2.4-3), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1max[0]) , '.2f'), format((p3max[0]) , '.2f'), format((p2max[0]) , '.2f')), fontsize=17, color='red')
            
            plt.text(0.0003, 10**(2.1-3), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1mid[0]) , '.2f'), format((p3mid[0]) , '.2f'), format((p2mid[0]) , '.2f')), fontsize=17, color='green')
                        
            plt.text(0.0003, 10**(1.8-3), r'$\rho \propto R^{%s, %s, %s}$' %(format((p1min[0]) , '.2f'), format((p3min[0]) , '.2f'), format((p2min[0]) , '.2f')), fontsize=17, color='blue')

            txt="Caption: Density profile of all the 41 globular clusters in the last snapshot of the simulation. The black distribution correspond to the \n density profile for all the  particles that belong to Globular Clusters. The blue, green and red represent the \n most negative energies, mid energies and energies in the peak of the distribution,  respectively. \n The lines correspond to the linear fit to the double power law in each case. The convertion from parsecs to degrees was implemented assuming \n the Galactocentric distance R0 = 8.178 ± 13stat. ± 22sys. pc.  (The GRAVITY Collaboration, A&A 625, L10 (2019))."
            plt.text(0.2, 0.0004, txt, wrap=True, horizontalalignment='center', fontsize=12)
            
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
