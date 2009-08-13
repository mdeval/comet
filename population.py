#! /usr/bin/python
"""
Plot level population as a function of radius
Reads output data from amc
"""

import sys
import matplotlib.pyplot as plt
import numpy
import fileinput 
import re 
import comet

# find number of lines to skip
for line in fileinput.input(sys.argv[1]):
	if re.match("@",line):
		skip=fileinput.lineno()
		break
fileinput.close()

r = numpy.loadtxt(sys.argv[1],skiprows=skip,comments="#",usecols=(1,))
rmax = numpy.loadtxt(sys.argv[1],skiprows=skip,comments="#",usecols=(2,))
lp = numpy.loadtxt(sys.argv[1],skiprows=skip,usecols=(range(9,9+len(comet.levels))))

for i in range(len(comet.levels)): # plot level population
	plt.loglog(numpy.append(r,rmax[-1]),numpy.append(lp[:,i],lp[-1,i]))
	plt.text(r[5],lp[5,i],"%d$_{%d%d}$" % (comet.levels[i][0],comet.levels[i][1],comet.levels[i][2]))
# 	plt.text(1.1*r[-1],lp[-1,i],"%d$_{%d%d}$" % (comet.levels[i][0],comet.levels[i][1],comet.levels[i][2]),\
# 			ha='left',va='center')
plt.xlabel('r [m]')
plt.ylabel('relative population')
# plt.axis([1e4,1e9,1e-5,1])
plt.show()
