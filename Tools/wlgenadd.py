#!/usr/bin/python2
#This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
#Copyright (c) 2023 Gururaj A. Wagle.  All rights reserved.
"""
Created on Tue Sep 13 10:33:21 2022

@author: Gururaj Wagle

Create a input.wlgrid file for a custom grid. Same can be used to create
input.fluxwl file (by modifying open statement at line 15).
"""
import numpy as np
from matplotlib import pyplot as plt

# open input.wlgrid file
f = open('input.wlgrid','w')
print("File name: ",f.name)

# construct grid
## set min, max values and number of bins
ng = 50
ngsub = 2000
wlmin = 1e-8
wlmax = 50e-5
wlgrid = np.zeros(ng)
for i in np.arange(ng):
    wlgrid[i] = (wlmin*(wlmax/wlmin)**(float(i)/float(ng)))
#print(wlgrid)

fig99 = plt.figure(num=99)
fig99.clear()
ax99 = fig99.add_subplot(1,1,1)
ax99.plot(wlgrid[:-1],np.diff(wlgrid))

wlgrid1 = wlgrid[wlgrid<3e-5]
wlgrid2 = wlgrid[wlgrid>10e-5]

wlgrid_res = np.logspace(np.log10(3e-5),np.log10(10e-5),1500)

wlgrid = np.concatenate((wlgrid1,wlgrid_res,wlgrid2))

ax99.plot(wlgrid[:-1],np.diff(wlgrid))
fig99.tight_layout()

#print(wlgrid)
# convert grid values to string
f.write(str(len(wlgrid))+" ")
for i in range(0,len(wlgrid)):
    f.write("%.6e "%wlgrid[i])
f.writelines("#")

# close input.wlgrid file
f.close()
