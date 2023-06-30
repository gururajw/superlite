#!/usr/bin/env python

#-- load modules
import time
import h5py
import numpy as np
import argparse
import matplotlib.pyplot as plt

#-- constructing parse field
parser = argparse.ArgumentParser(description = 'Coarsen HDF5 supernova data and convert to ASCII.')
parser.add_argument('f', type=str, help='HDF5 file.')

#-- parsing arguments
args = parser.parse_args()

#-- make an element dictionary of el/isotope to column index
el_dict = {'h1  ':0, 'he3 ':1, 'he4 ':1, 'c12 ':2, 'n14 ':3, 'o16 ':4,
           'ne20':5, 'mg24':6, 'si28':7, 's32 ':8, 'ar36':9,
           'ca40':10, 'ti44':11, 'cr48':12, 'fe52':13, 'fe54':13,
           'ni56':14}
inv_el_dict = {0:'H', 1:'He', 2:'C', 3:'N', 4:'O',
               5:'Ne', 6:'Mg', 7:'Si', 8:'S', 9:'Ar',
               10:'Ca', 11:'Ti', 12:'Cr', 13:'Fe', 14:'Ni'}
nelem = len(inv_el_dict)

#-- parse file
f = h5py.File(args.f, 'r')

#-- print out keys and group properties
keys = list(f.keys())
print()
print(len(keys), ' keys in file:')
for key in keys:
    print()
    print(' ', key, ', shape = ', f[key].shape)

#-- get the number of blocks
nblock = np.size(f['node type'])
print()
print('nblock = ', nblock)

#-- try averaging over blocks
nx_z = 16
ny_z = 16
nz_z = 16
nx_z_inv = 1.0 / nx_z
ny_z_inv = 1.0 / ny_z
nz_z_inv = 1.0 / nz_z

#-- time coarsening step
start = time.time()

coor = []
vels = []
mass = []
temp = []
dens = []
abun = []
for iblock in range(nblock):
    print(iblock, ' / ', nblock)
    if(f['node type'][iblock] == 1):
        #-- block edge coordinates
        x_min = f['bounding box'][iblock, 0, 0]
        x_max = f['bounding box'][iblock, 0, 1]
        if (x_max < x_min): print('x_max < x_min for iblock = ', iblock)
        y_min = f['bounding box'][iblock, 1, 0]
        y_max = f['bounding box'][iblock, 1, 1]
        if (y_max < y_min): print('y_max < y_min for iblock = ', iblock)
        z_min = f['bounding box'][iblock, 2, 0]
        z_max = f['bounding box'][iblock, 2, 1]
        if (z_max < z_min): print('z_max < z_min for iblock = ', iblock)
        coor.append(np.array([x_min, x_max, y_min, y_max, z_min, z_max]))
        #-- block edge velocities
        velx_min = np.mean(f['velx'][iblock, 0, :, :])
        velx_max = np.mean(f['velx'][iblock, nx_z - 1, :, :])
        vely_min = np.mean(f['vely'][iblock, :, 0, :])
        vely_max = np.mean(f['vely'][iblock, :, ny_z - 1, :])
        velz_min = np.mean(f['velz'][iblock, :, :, 0])
        velz_max = np.mean(f['velz'][iblock, :, :, nz_z - 1])
        vels.append(np.array([velx_min, velx_max, vely_min, vely_max, velz_min, velz_max]))
        #-- total block mass
        dx = nx_z_inv * (x_max - x_min)
        dy = ny_z_inv * (y_max - y_min)
        dz = nz_z_inv * (z_max - z_min)
        mass_ib = dx * dy * dz * np.sum(f['dens'][iblock, :, :, :])
        mass.append(mass_ib)
        dens.append(mass_ib / ((x_max - x_min) * (y_max - y_min) * (z_max - z_min)))
        #-- volume average temperature
        temp.append(np.mean(f['temp'][iblock, :, :, :]))
        #-- mass average abundances
        mfrac = np.array(nelem * [0.0])
        for el_key, el_val in el_dict.items():
            mfrac[el_val] += dx * dy * dz * np.sum(f['dens'][iblock, :, :, :] * f[el_key][iblock, :, :, :])
        mfrac /= mass_ib
        if (np.sum(mfrac) != 1.0): print('sum(mfrac) = ', np.sum(mfrac), ' for iblock = ', iblock)
        mfrac[mfrac < 1e-15] = 0.0
        abun.append(mfrac)

#-- time coarsening step
end = time.time()

print('cpu time = ', end - start)

coor = np.array(coor)
vels = np.array(vels)
mass = np.array(mass)
temp = np.array(temp)
dens = np.array(dens)
abun = np.array(abun)
print(np.sum(mass) / 2e33, ' Msol')

#-- plot
radi = np.array([np.linalg.norm(cr) for cr in coor])
velr = np.array([np.linalg.norm(vr) for vr in vels])

plt.plot(radi, dens, '*')
plt.xlabel('Radius [cm]')
plt.ylabel('Density [g/cm$^3$]')
plt.savefig('rho_hdf5.png')
plt.clf()
#plt.show()

plt.plot(radi, velr, '*')
plt.xlabel('Radius [cm]')
plt.ylabel('Speed [cm/s]')
plt.savefig('speed_hdf5.png')
plt.clf()
#plt.show()

plt.plot(radi, temp, '*')
plt.xlabel('Radius [cm]')
plt.ylabel('Temperature [K]')
plt.savefig('temp_hdf5.png')
plt.clf()
#plt.show()

#-- save the processed data
raw_str = np.array(np.size(mass) * [(15 + nelem) * [0.0]])
raw_str[:, :6] = coor
raw_str[:, 6:12] = vels
raw_str[:, 12] = mass
raw_str[:, 13] = temp
raw_str[:, 14] = dens
raw_str[:, 15:] = abun
#-- generate file headers
nx = np.size(np.unique(raw_str[:,0]))
if (nx != np.size(np.unique(raw_str[:,1]))): print('Invalid nx')
ny = np.size(np.unique(raw_str[:,2]))
if (ny != np.size(np.unique(raw_str[:,3]))): print('Invalid ny')
nz = np.size(np.unique(raw_str[:,4]))
if (ny != np.size(np.unique(raw_str[:,5]))): print('Invalid nz')
hd1 = 'cartesian\n'
hd2 = '  '+str(nx)+'     '+str(ny)+'     '+str(nz)+'     '+str(15 + nelem)+'     '+str(nelem)+'\n'
hd3 = '   x_left'
expllabs=['x_right','y_left','y_right','z_left','z_right',
          'vx_left', 'vx_right','vy_left','vy_right','vz_left',
          'vz_right','mass','temp','dens']
for lbl in expllabs: hd3+=lbl.rjust(12)
for el_key, el_val in inv_el_dict.items(): hd3+=el_val.rjust(12)
hd = hd1 + hd2 + hd3
#-- save
np.savetxt('input.allstr_x'+str(nx)+'y'+str(ny)+'z'+str(nz), raw_str, fmt='% .4e',delimiter=' ', header=hd)
