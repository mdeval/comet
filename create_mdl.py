#! /usr/bin/python
"""
Creates comet.mdl input data file

columns=id,ra,rb,ne,nm,tk,te,db,vr
id, inner radius, outer radius, electron number density, molecule number density,
kinematic temperature, electron temperature, Doppler b-coefficient, radial velocity
"""

import math
import comet
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--nshells", dest="nshells", default=50, type=int, help="number of radial shells")
(options, args) = parser.parse_args()

QH2O = 1e30 # s-1
vexp = 0.85 # km s-1
beta = 1.042e-5 # s-1
xre = 1.
xne = 1.
Rcs = 1125.*xre*(QH2O*1e-29)**.75 # km
Rrec = 3200.*xre*math.sqrt(QH2O*1e-29) # km
rh = 1. # AU
kion = 4.1e-7 # s-1

print """rmax=1E9
ncell=%d
tcmb=2.728
columns=id,ra,rb,nh,ne,nm,tk,te,db,vr
gas:dust=100
@""" % options.nshells

ra, rb = 0., 10.
for i in range(options.nshells):
    ra , rb = rb, rb*10**(5./options.nshells)
    r = (ra + rb)/2
    nH2O = QH2O /(4*math.pi*r**2*vexp)*math.exp(-r*beta/vexp)*.75*1e-15
    if r < Rcs:
        Te = comet.Tkin
    elif r > 2*Rcs:
        Te = 1e4
    else:
        Te = comet.Tkin + (1e4 - comet.Tkin)*(r/Rcs-1)
    krec = math.sqrt(300/Te)*7e-22 # s-1 km3
    ne = xne*math.sqrt(QH2O*kion/vexp/krec)/rh*(Te/300)**.15*(Rrec/r**2)*(1.-math.exp(-r/Rrec))*1e-15 + 5
    print i+1, ra*1e3, rb*1e3, nH2O, ne, nH2O, comet.Tkin, Te, 0.12012, vexp
