#! /usr/bin/python
"""
Plot H2O, electron densities and gas temperature.
Uses pyplot.twinx to have two different y scales.
"""

import numpy
import math
import matplotlib.pyplot as plt

rexp = numpy.arange(1,6,0.01)
r =10**rexp

QH2O = 1e29 # s-1
vexp = 0.85 # km/s
beta = 1.042e-5 # s-1
Tkin = 50
xre = 1
xne = 1
Rcs = 1125.*xre*(QH2O*1e-29)**.75 # km
Rrec = 3200.*xre*math.sqrt(QH2O*1e-29) # km
#rh = 1.49597870662e8 # AU in km
rh = 1 # AU
kion = 4.1e-7 # s-1
condlist = [ r < Rcs ]
Te = numpy.select (condlist, [Tkin])
condlist = [ r > 2*Rcs ]
Te = Te + numpy.select (condlist, [1e4])
condlist = [ Te < 30. ]
Te = Te + numpy.select (condlist,[Tkin + (1e4 - Tkin)*(r/Rcs-1)])

nH2O = QH2O /(4*math.pi*r**2*vexp)*numpy.exp(-r*beta/vexp)
krec = numpy.sqrt(300/Te)*7e-22 # s-1km3
ne = xne*numpy.sqrt(QH2O*kion/vexp/krec)/rh*(Te/300)**.15*(Rrec/r**2)*(1.-numpy.exp(-r/Rrec))
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.loglog(r,.75*nH2O*1e-15)
ax1.loglog(r,ne*1e-15+5,'c-.') # include electrons from the solar wind
ax1.set_xlabel("r [km]")
ax1.set_ylabel(r'n [cm$^{-3}$]')
ax1.set_ylim([1,1e12])

ax2 = ax1.twinx()
ax2.loglog(r,Te,'r--')
ax2.loglog([10,1e6],[Tkin,Tkin],'g:')
ax2.set_ylim([1,1e5])
ax2.set_ylabel("Te [K]")
plt.show()
