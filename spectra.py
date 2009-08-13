#! /usr/bin/python
"""
Plot line profiles.

Reads output data file from miriad. 
"""

import sys
import os
import re
from optparse import OptionParser
import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
import comet

def create_sky(filename, transtr, pixel):
    """creates sky.inp"""
    f = open(filename+".inp",'w')
    f.write( """source=%s.pop
    outfile=%s
    trans=%s
    pix=%d,%.02f
    chan=51,0.1
    distance=4.84825e-6
    units=K
    go
    q\n""" % ((filename,)*2+(transtr[0:-1],int(maxres/pixel+1),pixel)))
    f.close()

parser = OptionParser()
parser.add_option("-c", "--convol", action="store_true", dest="convol",
        default=False, help="run Miriad's convol")
parser.add_option("-f", "--file", dest="file", help="file")
(options, args) = parser.parse_args()

rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 8
# rcParams['text.fontsize'] = 6        #deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6

trans = (1,)
maxres=2.5*max([comet.HFWHM[i] for i in list(trans)])
transtr = ""
pixel = 0.05 # arcsec pixel size
for i in trans: transtr+=str(i+1)+","
# Creat sky.inp
create_sky(options.file, transtr, pixel)
# Run sky tracer program
os.system("rm -rf %s_00?" % options.file)
os.system("sky %s.inp" % options.file)

skip=4
exp=re.search("1e(\d\d)",options.file).group(1)
for i in trans: # loop over transitions
    up = comet.transitions[i][0]
    low = comet.transitions[i][1]
    boxsize = 1.25*comet.HFWHM[i] # arcsec
    filename = options.file+"_%03d"%int(i+1)
    print filename
    scale = pixel**2/comet.HFWHM[i]**2/1.13309
    if options.convol:
        os.system("rm -rf %s.cnv"%filename)
        os.system("convol map=%s out=%s.cnv fwhm=%f scale=%f"%(filename,filename,comet.HFWHM[i],scale))
        os.system("imspect in=%s.cnv region='arcsec,box(%d,-%d,-%d,%d)' log=%s.txt"%
                ((filename,)+(boxsize,)*4+(filename,)) )
    else:
        os.system("imspect in=%s region='arcsec,box(%d,-%d,-%d,%d)' log=%s.txt"%
                ((filename,)+(boxsize,)*4+(filename,)) )
    v = numpy.loadtxt(filename+".txt",skiprows=skip,usecols=(1,))
    t_B = numpy.loadtxt(filename+".txt",skiprows=skip,usecols=(2,))

    plt.figure(figsize=(3.5,3))
    plt.plot(v,t_B)
    plt.xlabel('v [km s$^{-1}$]')
    plt.ylabel('average intensity [K]')
    plt.xlim([v[0],v[-1]])
#     plt.xlim([-2.5,2.5])
    plt.text(v[int(v.size/5.)],3*plt.ylim()[1]/4.,"%d$_{%d%d}$-%d$_{%d%d}$" % (comet.levels[up][0],comet.levels[up][1],
        comet.levels[up][2],comet.levels[low][0],comet.levels[low][1],comet.levels[low][2]),fontsize='large')
    plt.text(v[int(v.size/10.)],2.5*plt.ylim()[1]/4.,"Q$_{H_2O}$ = 10$^{%s}$ s$^{-1}$" % exp,fontsize='large')
    plt.savefig(filename+".png")
    plt.close()
