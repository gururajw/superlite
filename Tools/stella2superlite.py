#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 11:52:10 2021

@author: Gururaj Wagle

An input.str file for SuperLite Snapshot (1d) is created from Stella profile
file

Example:
>> python stella_to_supernu_lite_IIn.py -d=20 -t -s --show-plots
"""

import os,re,sys,shutil
import argparse
import numpy as np
from myfuncts import hasNum,find_nearest_ind,safe_log_array,replace_line
from myfuncts import load_stella_prof
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl

home = os.path.expanduser('~')
M_sun = 1.98841e33 #grams

'''
1a. Input parameters
'''
parser = argparse.ArgumentParser(
    description='Convert Stella profile file to SuperLite input.str',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o','--out-path',type=str,default=os.getcwd(),
                    help='Absolute path to the folder to save input.str'
                    ' file')
parser.add_argument('-p','--path',type=str,default=os.getcwd(),
                    help='Absolute path to the Stella folder'
                    ' containing mesa.dayNNN_post_Lbol_max.data file')
parser.add_argument('-d','--day',type=float,default=30,
                    help='The day since L_{bol,max} to create a snapshot')
parser.add_argument('-t','--truncate',action='store_true',
                    help='Truncate the profile for tau>tau_thresh')
parser.add_argument('--tau-thresh',type=float,default=100,
                    help='The Stella calculated optical depth at which'
                    'the profile is truncated')
parser.add_argument('-s','--sanity-check',action='store_true',
                    help='Check for consistency')
parser.add_argument('-n','--renorm',action='store_true',
                    help='Whether to renormalize the mass fractions in each'
                    'cell if they do not add up to 1')
parser.add_argument('--homologous-vel',action='store_true',
                    help='Create homologize the velocity profile')
parser.add_argument('--inter',action='store_true',
                    help='Interpolate the original data to a higher'
                    ' resolution for 10> tau > 0.1')
parser.add_argument('--merge',action='store_true',
                    help='Merge "n-merge" cells interior to tau=10 into 1 cell')
parser.add_argument('--n-merge',type=int,default=5,
                    help='Number of inner cells to merge into 1')
parser.add_argument('-i','--only-info',action='store_true',
                    help='Print only the model info without creating str file,'
                    'and display the profile plots')
parser.add_argument('-v','--verbose',action='store_true',
                    help='increase terminal output verbosity')
parser.add_argument('--show-plots',action='store_true',
                    help='display the profile plots before saving')
parser.add_argument('--out-prefix',type=str,default='Model',
                    help='Common prefix for output profile plots')

### Esssential inputs
args = parser.parse_args()
stella_profile_path = args.path
outdir = args.out_path
day = args.day
out_prefix = args.out_prefix

if day>=1:
    day = int(day)
if day<10:
    day = '00'+str(day)
elif 10<=day<100:
    day = '0'+str(day)
elif day>=100:
    day = str(day)

### Set flags
sanity_check = args.sanity_check
renormalize= args.renorm
homologous_vel = args.homologous_vel
truncate = args.truncate
interpolate = args.inter
merge_inner = args.merge
only_info = args.only_info
verbose = args.verbose
show = args.show_plots

### set values
tau_thresh= args.tau_thresh
nmerge = args.merge

### Plot settings
mpl.rc('figure',autolayout=True)
mpl.rc('xtick', labelsize=14, top=True, direction='inout')
mpl.rc('xtick.major',size=12)
mpl.rc('xtick.minor',visible=True, size=8)
mpl.rc('ytick',labelsize=14,right=True,direction='inout')
mpl.rc('ytick.major',size=12)
mpl.rc('ytick.minor',visible=True, size=8)
mpl.rc('font',size=14)
mpl.rc('axes',labelsize=14,linewidth=0.8)
mpl.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.cividis.colors)
mpl.rc('legend',frameon=False)
mpl.rc('lines',lw=2)

## Colorblind-friendly colors
plt.style.use('tableau-colorblind10')
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors = ['#1F77B4','#FF7F0E','#2CA02C','#D62728','#9467BD',
          '#8C564B','#E377C2','#7F7F7F','#BCBD22','#17BECF']

"""
1b. Read Stella Profile file
"""
stella_file = stella_profile_path+'/mesa.day'+day+'_post_Lbol_max.data'

if not os.path.exists(stella_file):
    print("The required Stella profile file doesn't exist:\n%s\n"
          "Use the option -p/--path to set correct path or -d/--day to "
          "set correct day for which the profile file exists." % stella_file)
    sys.exit()

if not os.path.exists(outdir):
    mkdir = input("The directory specified for saving the input.str file"
                  " doesn't exist.\nCreate a new directory at %s?\nEnter Yes "
                  "to continue, anything else to quit: " % outdir)
    if mkdir.lower()=='yes':
        os.mkdir(outdir)
    else:
        sys.exit()

data_stella = load_stella_prof(stella_file)

if not (verbose or only_info):
    f = open(outdir+'/'+out_prefix+'_profile_day'+day+'.info',"w")
    temp = sys.stdout
    sys.stdout = f

L_bol =  data_stella['Lum'][-1]
R_tau_23 = data_stella['r_center_cm'][
    find_nearest_ind(data_stella['tau'],2/3)]


print("--------------------------------------------------------")
print("TOTAL MASS IN THE INPUT MODEL")
print("--------------------------------------------------------")
print("Total Nickel mass in Stella model: %.4e Msun" %
  ((data_stella['ni56']*data_stella['m_cell_g']).sum()/M_sun))
print("Ejecta Mass in Stella model: %.2f Msun\n"
      % (data_stella['m_cell_g'].sum()/M_sun))
print("--------------------------------------------------------")
print("Cell center radii, inner and outer: %.4e, %.4e cm" %
      (data_stella['r_center_cm'][0],data_stella['r_outer_cm'][-1]))
print("L_{bol}: %.4e erg/s\n" % L_bol)
print("R_{tau=2/3} is at %.4e cm" % R_tau_23)
print("R_{tau=1} is at %.4e cm"
      % (data_stella['r_center_cm'][
          find_nearest_ind(data_stella['tau'],1)]))
print("R_{tau=3} is at %.4e cm"
      % (data_stella['r_center_cm'][
          find_nearest_ind(data_stella['tau'],3)]))
print("--------------------------------------------------------")

f_vel = interp1d(data_stella['r_center_cm'],data_stella['v_center_cmps'],
                 kind = 'slinear',fill_value="extrapolate")
vel_right = f_vel(data_stella['r_outer_cm'])
### Plot several profiles for sanity check
ind_tau1 = find_nearest_ind(data_stella['tau'],1)

fig1 = plt.figure(figsize=(16,12))
fig1.clear()
plt.suptitle('day = %s' % day)
ax11 = fig1.add_subplot(2,2,1)
ax11.set_title('Velocity')
ax11.set_xlabel('$ \\rm R_{outer} \ (\\times 10^{14} \ cm)$')
ax11.set_ylabel('$ \\rm v_{outer} \ (\\times 10^3 \ km \ s^{-1})$')
ax11.plot(data_stella['r_outer_cm']/1e14,vel_right/1e8,
          color=colors[0],label="Stella")
ax11.axvline(data_stella['r_center_cm'][ind_tau1]/1e14,
             color='darkgrey',ls='--')
ax11.text(data_stella['r_center_cm'][ind_tau1]/1e14,
          ax11.get_ylim()[1]*.5,'$\\rm \\tau_{Stella} = 1$',
          ha='right',rotation=90,fontsize=12)
ax11.legend()
ax12 = fig1.add_subplot(2,2,2)
ax12.set_title('Cumulative Mass')
ax12.set_xlabel('$ \\rm R_{outer} \ (\\times 10^{14} \ cm)$')
ax12.set_ylabel('Mass ($\\rm M_{\\odot}$)')
ax12.plot(data_stella['r_outer_cm'][1:]/1e14,
          np.cumsum(data_stella['m_cell_g'][1:]/M_sun),
          color=colors[0],label="Stella")
ax12.axvline(data_stella['r_center_cm'][ind_tau1]/1e14,
             color='darkgrey',ls='--')
ax12.text(data_stella['r_center_cm'][ind_tau1]/1e14,
          ax12.get_ylim()[1]*.5,'$\\rm \\tau_{Stella} = 1$',
          ha='right',rotation=90,fontsize=12)
ax12.legend()
ax13 = fig1.add_subplot(2,2,3)
ax13.set_title('Temperature')
ax13.set_xlabel('$ \\rm R_{outer} \ (\\times 10^{14} \ cm)$')
ax13.set_ylabel('$\\rm T \ (\\times 10^4 \ K)$')
ax13.plot(data_stella['r_outer_cm']/1e14,data_stella['Tavg']/1e4,
          color=colors[0],label="Stella,T$_{avg}$")
ax13.plot(data_stella['r_outer_cm']/1e14,data_stella['Trad']/1e4,
          color=colors[1],label="Stella,T$_{rad}$")
ax13.axvline(data_stella['r_center_cm'][ind_tau1]/1e14,
             color='darkgrey',ls='--')
ax13.text(data_stella['r_center_cm'][ind_tau1]/1e14,
          ax13.get_ylim()[1]*.5,'$\\rm \\tau_{Stella} = 1$',
          ha='right',rotation=90,fontsize=12)
ax13.legend()
ax14 = fig1.add_subplot(2,2,4)
ax14.set_title('Optical Depth')
ax14.set_xlabel('$ \\rm R_{outer} \ (\\times 10^{14} \ cm)$')
ax14.set_ylabel('$\\rm log_{10}\\tau$')
ax14.plot(data_stella['r_outer_cm']/1e14,safe_log_array(data_stella['tau']),
          color=colors[0],label="Stella")
ax14.axvline(data_stella['r_center_cm'][ind_tau1]/1e14,
             color='darkgrey',ls='--')
ax14.text(data_stella['r_center_cm'][ind_tau1]/1e14,
          ax14.get_ylim()[1]*.5,'$\\rm \\tau_{Stella} = 1$',
          ha='right',rotation=90,fontsize=12)
ax14.legend()

if outdir==os.getcwd():
    print("\nWARNING - Output directory is not specified. The input.str file"
          " will be created at the current location. This will overwrite any"
          " existing input.str file")

data_stella_original = data_stella.copy() # back up original Stella data

if truncate:
    ind_tau_thresh = find_nearest_ind(data_stella['tau'],tau_thresh)
    data_stella = data_stella[ind_tau_thresh:]
    if verbose:
        print("\nTruncating the profile at tau=100, R=%.4e\n" %
              data_stella['r_center_cm'][0])
else:
    ind_tau_thresh = -1

"""
2. Add isotopes of similar type into elemental mass fractions from profile 0
"""
### Elements data type (as input to SuperLite)
dt_elem = [('h',float),('he',float),('c',float),('n',float),('o',float),
           ('ne',float),('na',float),('mg',float),('si',float),('s',float),
           ('ar',float),('ca',float),('ti',float),('cr',float),('fe',float),
           ('co',float),('ni',float),('cr48',float),('fe52',float),
           ('co56',float),('ni56',float)]

### Array for elemental mass fractions
abund=np.zeros(data_stella.shape,dtype=dt_elem)

# save isotopic abundances for radioactive isotopes cr48, fe52, co56, & ni56
for name in abund.dtype.names:
    if (hasNum(name) and name in data_stella.dtype.names):
        abund[name] = data_stella[name]

### Add isotopic abundances for each element to elemental abundances
niso=0
for name in data_stella.dtype.names:
    if hasNum(name):
        niso+=1
        (elem,iso)=([re.findall(r'(\w+?)(\d+)',name)[0]][0])
        if elem in abund.dtype.names:
            abund[elem]+=data_stella[name]

### Sanity check: sum of mass fractions = 1 in all cells
sum_mfrac=np.zeros(len(abund))
for name in abund.dtype.names:
    if not(name =='ni56' or name=='co56' or name=='fe52' or name =='cr48'):
        for i in range(0,len(abund)):
            sum_mfrac[i]+=abund[name][i]

if sanity_check:
    print("--------------------------------------------------------")
    print("CHECKING NORMALIZATION OF MASS FRACTIONS IN THE INPUT")
    print("--------------------------------------------------------")
    print("The mean total mass fraction per zone is:%.5f" % np.mean(sum_mfrac))
    print("The median total mass fraction per zone is:%.5f"
          % np.median(sum_mfrac))
    print("The mean excess/deficit from 1.0 is:%.2e\n"
          % np.mean(sum_mfrac - 1.0))

    if np.mean(sum_mfrac - 1.0) > 1.e-4:
        print("The final mass fractions do not normalize to unity"
               " within 10^-4.")
        if renormalize:
            print("renormalizing to unity...\n")


### Renormalize the abundances, if the flag is set to True
if (renormalize and np.mean(sum_mfrac - 1.0) > 1e-4):
    for name in abund.dtype.names:
        for i in range(0,len(abund)):
            abund[name][i]/=sum_mfrac[i]

### Sanity check: post-renormalization sum of mass fractions = 1 in all cells
if (sanity_check and renormalize and np.mean(sum_mfrac - 1.0) > 1e-4):
    sum_mfrac=np.zeros(len(abund))
    for name in abund.dtype.names:
         if not(name =='ni56' or name=='co56' or name=='fe52' or name =='cr48'):
             for i in range(0,len(abund)):
                 sum_mfrac[i]+=abund[name][i]
    if verbose:
        print("After renormalization:")
        print("The mean total mass fraction per zone is:%.5f" % np.mean(sum_mfrac))
        print("The median total mass fraction per zone is:%.5f"
              % np.median(sum_mfrac))
        print("The mean excess/deficit from 1.0 is:%.2e\n"
              % np.mean(sum_mfrac - 1.0))

        if np.mean(sum_mfrac - 1.0) > 1.e-4:
            print("The final mass fractions do not normalize to unity"
                  " within 10^-4.\n")
        else:
            print("The final mass fractions renormalized to unity.\n")

if(sanity_check):
    print("--------------------------------------------------------")

### Create output data array
vel_right = np.zeros(len(data_stella['ind']))
if homologous_vel:
    print("\nCreating homologous velocity profile\n")
### Linear regression to velocity
#    slope,intercept = np.polyfit(data_stella['r_center_cm'],
#                                 data_stella['v_center_cmps'],1)
    ### intercept is not added (see Ryan's script in src/tools)
#    vel_new = slope*data_stella['r_center_cm']+intercept

    slope1 = ((data_stella['v_center_cmps'][-1] -
              data_stella['v_center_cmps'][0])/
              (data_stella['r_center_cm'][-1] -
               data_stella['r_center_cm'][0]))
    vel_new = slope1*data_stella['r_center_cm']
    intercept1 = vel_new[-1]- data_stella['v_center_cmps'][-1]
    vel_new = vel_new-intercept1
    for i in range(0,len(data_stella['ind'])):
        vel_right[i] = (data_stella['r_outer_cm'][i]*
                        vel_new[i]/
                        data_stella['r_center_cm'][i])
else:
    f_vel = interp1d(data_stella['r_center_cm'],data_stella['v_center_cmps'],
                     kind = 'slinear',fill_value="extrapolate")
    vel_right = f_vel(data_stella['r_outer_cm'])
    # for i in range(0,len(data_stella['ind'])):
    #     vel_right[i] = (data_stella['r_outer_cm'][i]*
    #                     data_stella['v_center_cmps'][i]/
    #                     data_stella['r_center_cm'][i])

n_zones = len(data_stella['ind'])-1

if interpolate:
    ind_tau10 = find_nearest_ind(data_stella['tau'],10)
    ind_tau0p1 = find_nearest_ind(data_stella['tau'],0.1)
    if (merge_inner and not truncate):
        n_merge = int(ind_tau10*nmerge) # New number of merged cells
        n_remain = ind_tau10 - n_merge*nmerge # Remainder cells
    else: # Don't merge
        n_merge = 0
        n_remain = ind_tau10

    n_interp = (ind_tau0p1 - ind_tau10)*nmerge
    n_outer = n_zones - ind_tau0p1
    n_zones = n_merge + n_remain + n_interp + n_outer
    if verbose:
        print("Interpolating data between tau = 10 and 0.1")
        print("Number of zones for tau>10: %d"%(n_merge+n_remain))
        print("Number of zones for 10<tau<0.1: %d"%(n_interp))
        print("Number of zones for tau<0.1: %d"%n_outer)
        print("New number of zones: %d"%n_zones)

    # calculate cumulative masses
    m_cum = np.cumsum(data_stella['m_cell_g'])
    m_cum_elem = np.zeros(len(data_stella['ind']),dtype=dt_elem)
    for name in abund.dtype.names:
        m_cum_elem[name]=np.cumsum(np.multiply(data_stella['m_cell_g'],
                                               abund[name]))
    # interpolation functions
    f_mass = interp1d(data_stella['r_outer_cm'],m_cum,
                      kind = 'slinear',fill_value="extrapolate")
    f_radtemp = interp1d(data_stella['r_center_cm'],data_stella['Trad'],
                         kind = 'slinear',fill_value="extrapolate")
    f_temp = interp1d(data_stella['r_center_cm'],data_stella['Tavg'],
                      kind = 'slinear',fill_value="extrapolate")

    # merge inner cells
    if merge_inner:
        xleft_merge = np.zeros(n_merge)
        xcenter_merge = np.zeros(n_merge)
        xright_merge = np.zeros(n_merge)
        vx_merge = np.zeros(n_merge)
        mass_merge = np.zeros(n_merge)
        radtemp_merge = np.zeros(n_merge)
        temp_merge = np.zeros(n_merge)
        m_elem_merge = np.zeros(n_merge,dtype=dt_elem)
        abund_merge = np.zeros(n_merge,dtype=dt_elem)
        for ind in range(0,n_merge):
            xright_merge[ind] = data_stella['r_outer_cm'][((ind+1)*4-1)]
            if ind==0:
                xleft_merge[ind] = data_stella['r_center_cm'][0]
            else:
                xleft_merge[ind] = xright_merge[ind-1]
            xcenter_merge[ind] = 0.5*(xright_merge[ind]+xleft_merge[ind])
            vx_merge[ind] = vel_right[(ind+1)*4-1]
            temp_merge[ind] = f_temp(xcenter_merge[ind])
            radtemp_merge[ind] = f_radtemp(xcenter_merge[ind])
            if ind>0:
                mass_merge[ind] = m_cum[((ind+1)*4)-1]-m_cum[ind*4-1]
                for name in abund.dtype.names:
                    m_elem_merge[name][ind] = (m_cum_elem[name][((ind+1)*4)-1]-
                                               m_cum_elem[name][ind*4-1])
            else:
                mass_merge[ind] = m_cum[((ind+1)*4)-1]
                for name in abund.dtype.names:
                    m_elem_merge[name][ind] = m_cum_elem[name][((ind+1)*4)-1]
            for name in abund.dtype.names:
                abund_merge[name][ind]=m_elem_merge[name][ind]/mass_merge[ind]

    # keep remaining cells un_merged
    if n_remain>0:
        xleft_remain = np.zeros(n_remain)
        if merge_inner:
            xleft_remain = data_stella['r_outer_cm'][(n_merge*4-1):
                                                     (n_merge*4+n_remain-1)]
        else:
            xleft_remain[0] = data_stella['r_center_cm'][0]
            xleft_remain[1:] = data_stella['r_outer_cm'][0:(n_merge*4+n_remain-1)]
        xright_remain = data_stella['r_outer_cm'][n_merge*4:(n_merge*4+n_remain)]
        vx_remain = vel_right[n_merge*4:(n_merge*4+n_remain)]
        mass_remain = data_stella['m_cell_g'][n_merge*4:(n_merge*4+n_remain)]
        radtemp_remain = data_stella['Trad'][n_merge*4:(n_merge*4+n_remain)]
        temp_remain = data_stella['Tavg'][n_merge*4:(n_merge*4+n_remain)]
        abund_remain = abund[:][n_merge*4:(n_merge*4+n_remain)]

    # increase resolution near photosphere (10>tau>0.1)
    xleft_interp = np.zeros(n_interp)
    xcenter_interp = np.zeros(n_interp)
    xright_interp = np.zeros(n_interp)
    vx_interp = np.zeros(n_interp)
    mass_interp = np.zeros(n_interp)
    mcum_interp = np.zeros(n_interp)
    radtemp_interp = np.zeros(n_interp)
    temp_interp = np.zeros(n_interp)
    abund_interp = np.zeros((n_interp),dtype=dt_elem)

    dx = np.diff(data_stella['r_outer_cm'])
    # increase resolution
    j = ind_tau10
    for i in range(ind_tau10,ind_tau0p1):
        k = i-j
        # print(i,j,k)
        for x in range(0,4):
            # print(k+x)
            if (k==0 and x==0):
                xleft_interp[k+x]=data_stella['r_outer_cm'][i-1]
            else:
                xleft_interp[k+x]=xright_interp[k+x-1]
            # print(k+x)
            xright_interp[k+x] = (data_stella['r_outer_cm'][i-1]+
                                  (x+1)*dx[i-1]/4)
            xcenter_interp[k+x] = 0.5*(xleft_interp[k+x]+xright_interp[k+x])
            vx_interp[k+x] = ((vel_right[i-1]*(data_stella['r_outer_cm'][i]
                                               -xright_interp[k+x])+
                               vel_right[i]*(xright_interp[k+x]
                                             -data_stella['r_outer_cm'][i-1]))
                              /dx[i-1])
            for name in abund.dtype.names:
                abund_interp[name][k+x]=abund[name][i]
            mcum_interp[k+x] = f_mass(xright_interp[k+x])
            radtemp_interp[k+x] = f_radtemp(xcenter_interp[k+x])
            temp_interp[k+x] = f_temp(xcenter_interp[k+x])
        j-=3

    mass_interp[0] = mcum_interp[0] - m_cum[ind_tau10-1]
    mass_interp[1:] = np.diff(mcum_interp)


    # keep outer cells unchanged
    xleft_outer = data_stella['r_outer_cm'][ind_tau0p1-1:-1]
    xright_outer = data_stella['r_outer_cm'][ind_tau0p1:]
    vx_outer = vel_right[ind_tau0p1:]
    mass_outer = data_stella['m_cell_g'][ind_tau0p1:]
    radtemp_outer = data_stella['Trad'][ind_tau0p1:]
    temp_outer = data_stella['Tavg'][ind_tau0p1:]
    abund_outer = abund[:][ind_tau0p1:]

### Output files settings
dt_slite = [('x_left',float),('x_right',float),('vx_left',float),
          ('vx_right',float),('mass',float),('rad_temp',float),('temp',float),
          ('h',float),('he',float),('c',float),('n',float),('o',float),
          ('ne',float),('na',float),('mg',float),('si',float),('s',float),
          ('ar',float),('ca',float),('ti',float),('cr',float),('fe',float),
          ('co',float),('ni',float),('cr48',float),('fe52',float),
          ('co56',float),('ni56',float)]

outdata=np.zeros(n_zones,dtype=dt_slite)
nabund=0
if interpolate:
    if merge_inner:
        outdata['x_left'] = np.concatenate((xleft_merge,xleft_remain,
                                            xleft_interp,xleft_outer))
        outdata['x_right'] = np.concatenate((xright_merge,xright_remain,
                              xright_interp,xright_outer))
        outdata['vx_right'] = np.concatenate((vx_merge,vx_remain,
                                              vx_interp,vx_outer))
        outdata['mass'] = np.concatenate((mass_merge,mass_remain,
                                          mass_interp,mass_outer))
        outdata['radtemp'] = np.concatenate((radtemp_merge,radtemp_remain,
                                              radtemp_interp,radtemp_outer))
        outdata['temp'] = np.concatenate((temp_merge,temp_remain,
                                          temp_interp,temp_outer))
        for name in abund.dtype.names:
            if abund[name].any():
                outdata[name]=np.concatenate((abund_merge[name],
                                              abund_remain[name],
                                              abund_interp[name],
                                              abund_outer[name]))
                nabund+=1
    else:
        outdata['x_left'] = np.concatenate((xleft_remain,xleft_interp,
                                            xleft_outer))
        outdata['x_right'] = np.concatenate((xright_remain,xright_interp,
                                             xright_outer))
        outdata['vx_right'] = np.concatenate((vx_remain,vx_interp,vx_outer))
        outdata['mass'] = np.concatenate((mass_remain,mass_interp,mass_outer))
        outdata['radtemp'] = np.concatenate((radtemp_remain,radtemp_interp,
                                              radtemp_outer))
        outdata['temp'] = np.concatenate((temp_remain,temp_interp,temp_outer))
        for name in abund.dtype.names:
            if abund[name].any():
                outdata[name]=np.concatenate((abund_remain[name],
                                              abund_interp[name],
                                              abund_outer[name]))
                nabund+=1
else:
    outdata['x_left'] = data_stella['r_outer_cm'][:-1]
    outdata['x_right'] = data_stella['r_outer_cm'][1:]
    outdata['vx_left'] = vel_right[:-1]
    outdata['vx_right'] = vel_right[1:]
    outdata['mass'] = data_stella['m_cell_g'][1:]
    outdata['radtemp'] = data_stella['Trad'][1:]
    outdata['temp'] = data_stella['Tavg'][1:]
    for name in abund.dtype.names:
        if abund[name].any():
            outdata[name]=abund[name][1:]
            nabund+=1

### Sanity check
if any(outdata['mass']<=0):
    print('Error: Mass is <=0 in the output')
    print(outdata['mass'][outdata['mass']<=0])
    exit()
if sanity_check:
    print("--------------------------------------------------------")
    print("TOTAL MASS IN THE OUTPUT FILE")
    print("--------------------------------------------------------")
    print("Total Nickel mass in the output: %.2e Msun" %
          ((outdata['ni']*outdata['mass']).sum()/M_sun))
    print("Ejecta Mass in the output: %.2f Msun\n"
          % (outdata['mass'].sum()/M_sun))
    print("--------------------------------------------------------")
    print("Cell center radii, inner and outer: %.4e, %.4e cm" %
      (data_stella['r_center_cm'][0],data_stella['r_outer_cm'][-1]))
    print("L_{bol}: %.4e" %
          (data_stella['Lum'][-1]-data_stella['Lum'][0]))
    print("--------------------------------------------------------\n")
    print("--------------------------------------------------------")
    print("CHECKING NORMALIZATION OF MASS FRACTIONS IN THE OUTPUT")
    print("--------------------------------------------------------")
    print("The mean total mass fraction per zone is:%.5f" %
          np.mean(sum_mfrac))
    print("The median total mass fraction per zone is:%.5f"
          % np.median(sum_mfrac))
    print("The mean excess/deficit from 1.0 is:%.2e\n"
          % np.mean(sum_mfrac - 1.0))
    print("--------------------------------------------------------")

### Plot several profiles
ax11.plot(outdata['x_right']/1e14,outdata['vx_right']/1e8,
          'k--',lw=4,label="SuperLite input")
ax11.axvline(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
             color='darkgrey',ls='--')
ax11.text(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
          ax11.get_ylim()[1]*.1,'$\\rm \\tau_{Stella} = %d$'
          % int(data_stella_original['tau'][ind_tau_thresh]),
          ha='right',rotation=90,fontsize=12)
ax11.legend()
ax12.plot(outdata['x_right']/1e14,
          (np.sum(data_stella_original['m_cell_g'][:ind_tau_thresh+1])
          +np.cumsum(outdata['mass']))/M_sun,
          'k--',lw=4,label="SuperLite input")
ax12.axvline(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
             color='darkgrey',ls='--')
ax12.text(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
          ax12.get_ylim()[1]*.1,'$\\rm \\tau_{Stella} = %d$'
          % int(data_stella_original['tau'][ind_tau_thresh]),
          ha='right',rotation=90,fontsize=12)
ax12.legend()
ax13.plot(outdata['x_right']/1e14,outdata['temp']/1e4,
          'k--',lw=4,label="SuperLite input,T$_e$")
ax13.plot(outdata['x_right']/1e14,outdata['radtemp']/1e4,
         'k:',lw=4,label="SuperLite,T$_{rad}$")
ax13.axvline(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
             color='darkgrey',ls='--')
ax13.text(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
          ax13.get_ylim()[1]*.1,'$\\rm \\tau_{Stella} = %d$'
          % int(data_stella_original['tau'][ind_tau_thresh]),
          ha='right',rotation=90,fontsize=12)
ax13.legend()
ax14.axvline(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
             color='darkgrey',ls='--')
ax14.text(data_stella_original['r_center_cm'][ind_tau_thresh]/1e14,
          ax14.get_ylim()[1]*.05,'$\\rm \\tau_{Stella} = %d$'
          % int(data_stella_original['tau'][ind_tau_thresh]),
          ha='right',rotation=90,fontsize=12)
# ax14.plot(outdata['x_right'],outdata['radtemp']/1e4,
#           color=colors[1],ls=':',linewidth=4,label="Output Data")
# ax14.legend()


# if interpolate:
#     fig2 = plt.figure()
#     fig2.clear()
#     ax21 = fig2.add_subplot(2,1,1)
#     ax22 = fig2.add_subplot(2,1,2)
#     ax21.set_title("Zoning")
#     ax21.set_xlabel("$ \\rm R_{outer} \ (cm)$")
#     ax21.set_ylabel("dr (cm)")
#     ax21.plot(outdata['x_right'][:-1],np.diff(outdata['x_right']),
#               'k',linewidth=2)
#     ax22.set_xlabel("$ \\rm R_{outer} \ (cm)$")
#     ax22.set_ylabel("$ \\rm \\tau$")
#     ax22.plot(data_stella['r_outer_cm'],data_stella['tau'],'k',linewidth=2)
#     ax22.set_ylim(-0.1,20)
#     fig2.tight_layout()

if show:
    plt.show()

if not (verbose or only_info):
    f.close()
    sys.stdout = temp

if only_info:
    plt.show()
    sys.exit()
else:
    fig1.savefig(outdir+'/'+out_prefix+'_profile_day'+day+'.png')


"""
4. Prepare the output file input.str to be used as an input to Supernu
"""
nx = n_zones
ny = 1
nz = 1
ncol=0
for name in outdata.dtype.names:
    if outdata[name].any():
        ncol+=1

### input.str
### Header lines
outfile = open(outdir+'/input.str','w')
print("The new input.str file is created at %s/\n" % outdir)
outfile.writelines("# 1D spherical\n")
outfile.writelines("#      "+str(nx)+"          "+str(ny)+" "+str(nz)+" "
                   +str(ncol)+" "+str(nabund)+"\n")

outfile.write("#     x_left")

for name in outdata.dtype.names[1:]:
    if outdata[name].any():
        spaces = (14-len(name))*' '
        outfile.write("%s%s" % (spaces,name))
outfile.write("\n")
### data
for i in range(0,n_zones):
    for name in outdata.dtype.names:
        if outdata[name].any():
            outfile.write("%.6e  " % (outdata[name][i]))
    outfile.write("\n")
outfile.close()

### update parameters in input.par
shutil.copy2(home+'/workspaces/superlite/input.par.lte',
             outdir+'/input.par.lte')
shutil.copy2(home+'/workspaces/superlite/input.par.nlte',
             outdir+'/input.par.nlte')

file = outdir+'/input.par.lte'
replace_line(file, 15, " in_ndim = %d, 1, 1" % n_zones)
text = ("%.3e" % L_bol).replace('e','d')
replace_line(file, 36, " in_L_bol = %s" % text)
if float(day)>=1:
    replace_line(file, 57, "in_name = '%s, %dd, LTE'"
                 % (out_prefix,int(day)))
else:
    replace_line(file, 57, "in_name = '%s, %.1fd, LTE'"
                 % (out_prefix,float(day)))

file = outdir+'/input.par.nlte'
replace_line(file, 15, " in_ndim = %d, 1, 1" % n_zones)
text = ("%.3e" % L_bol).replace('e','d')
replace_line(file, 36, " in_L_bol = %s" % text)
text = ("%.3e" % data_stella['r_center_cm'][
          find_nearest_ind(data_stella['tau'],2/3)]).replace('e','d')
replace_line(file, 56, " in_R_phot = %s" % text)
if float(day)>=1:
    replace_line(file, 59, "in_name = '%s, %dd, NLTE'"
             % (out_prefix,int(day)))
else:
    replace_line(file, 59, "in_name = '%s, %.1fd, LTE'"
                 % (out_prefix,float(day)))
