#! /usr/bin/python
"""
Plot H2O, electron densities and gas temperature.

Uses pyplot.twinx to have two different y scales.
"""

import numpy
import sys
import matplotlib.pyplot as plt
import comet
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--save", action="store_true", default=False, dest="save", help="save png")
(options, args) = parser.parse_args()

skip = comet.skipno(sys.argv[1])

ra,rb,ne,nm,tkin,te = numpy.loadtxt(sys.argv[1],skiprows=skip,usecols=(1,2,4,5,6,7),unpack=True)
r = (ra+rb)/2.

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.loglog(r,nm)
ax1.loglog(r,ne,'c-.') # include electrons from the solar wind
ax1.set_xlabel("r [km]")
ax1.set_ylabel(r'n [cm$^{-3}$]')
ax1.set_ylim([1,1e12])
ax1.text(r[int(r.size/10.)],ne[int(r.size/10.)],"n$_e$")
ax1.text(r[int(r.size/10.)],nm[int(r.size/10.)],"n$_{H_2O}$")

ax2 = ax1.twinx()
ax2.loglog(r,te,'r--')
ax2.loglog(r,tkin,'g:')
# ax2.loglog([1e4,1e9],[comet.Tkin,comet.Tkin],'g:')
ax2.set_ylim([1,1e5])
ax2.set_ylabel("Te [K]")
ax2.text(r[int(8*r.size/10.)],te[int(8*r.size/10.)],"T$_e$")
ax2.text(r[int(8*r.size/10.)],tkin[int(8*r.size/10.)],"T$_k$")
if options.save:
    plt.savefig(sys.argv[1]+'.png')
else:
    plt.show()
