#! /usr/bin/python

import numpy
import math
import physcon as pc
from scipy import special
import comet

temp = numpy.arange(comet.Tkin,1e4,3e4)
temp = numpy.append(temp,1e4)

weight = [9.0, 9.0, 15.0, 15.0, 21.0, 21.0, 21.0]
transitions = [[1,0], [2,0], [3,1], [3,2], [4,2], [5,3], [5,4], [6,2], [6,5]]
# H2O-e collision rates (Xie & Mumma 1992)
linestrength = [1.5, 1.5, 1.5, .0,8333, 1.7441, 0.3003, 2.1809, 0.9225, 2.5434] # King et al 1947
#sigma = 8*math.pi/3/k**2*(2*jp+1)*d**2 * linestrength[transitions.index(i)]/(2*j+1)/(2*jp+1)*math.log((k+kp)/(k-kp))
#nu = sigma * numpy.sqrt(8*pc.k_B*numpy.array(temp)/math.pi/pc.m_e) # rate coefficients

# H2O-H2O collision rates for 100K
sigmah2o = numpy.array([2.924e-18,3.248e-18,2.162e-18,2.327e-18,2.627e-18,6.840e-19,2.204e-18,6.570e-19,2.322e-18]) # m2 (Zakharov et al. 2007)
Ch2o = 4.831e7*sigmah2o*math.sqrt(comet.TH2O) # [cm3 s-1] Buffa et al. 2000 (eq. 12)

# H2O-e collision rates (Zakharov et al. 2007) sec. 3.3
Aeins = numpy.array([ 3.458e-03, 5.593e-02, 2.564e-01, 3.058e-02, 5.048e-02, 2.634e-03, 1.648e-02, 3.318e-01, 2.288e-02]) # s-1
freq = 1e9*numpy.array([ 556.936, 1669.904, 2773.977, 1661.007, 1716.770, 1153.127, 1097.365, 3977.047, 1162.912 ]) # Hz
sigmaeij = 1e4*pc.m_e*pc.e**2*pc.c**3*Aeins/16/math.pi**2/pc.eps_0/pc.h**2/freq**4 # cm2
ve = 1e2*numpy.sqrt(8*pc.k_B*temp/math.pi/pc.m_e) # [cm s-1] averaged thermal speed

# approximate H2O-e collision rates value (Bensch & Bergin 2004)
Ce = 2.e-7 * numpy.ones(temp.size)

# create h2o.dat
print """!MOLECULE
Ortho H2O
!MOLECULAR WEIGHT
18.0
!NUMBER OF ENERGY LEVELS
7
!LEVEL + ENERGIES(cm^-1) + WEIGHT + J
  1    23.794356     9.0      1_0_1
  2    42.371741     9.0      1_1_0
  3    79.496382    15.0      2_1_2
  4   134.901638    15.0      2_2_1
  5   136.761650    21.0      3_0_3
  6   173.365803    21.0      3_1_2
  7   212.156362    21.0      3_2_1
!NUMBER OF RADIATIVE TRANSITIONS
9
!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)
    1     2     1  3.458e-03         556.936        26.7
    2     3     1  5.593e-02        1669.904        80.2
    3     4     2  2.564e-01        2773.977       159.9
    4     4     3  3.058e-02        1661.007       159.9
    5     5     3  5.048e-02        1716.770       162.5
    6     6     4  2.634e-03        1153.127       215.2
    7     6     5  1.648e-02        1097.365       215.2
    8     7     3  3.318e-01        3977.047       271.0
    9     7     6  2.288e-02        1162.912       271.0
!NUMBER OF COLL PARTNERS
2
! COLLISION PARTNER
2 oH2O - oH2O from Buffa et al. 2000
! NUMBER OF COLLISIONAL TRANSITIONS
9
! NUMBER OF COLLISION TEMPERATURES
1
! COLLISION TEMPERATURES
%.1f
! TRANS + UP + LOW + RATE COEFFS(cm^3 s^-1)""" % comet.Tkin
for i in transitions:
	print "%d %d %d %.3e" % (transitions.index(i)+1, i[0]+1, i[1]+1, Ch2o[transitions.index(i)])
print """!COLLISIONS BETWEEN  
3  H2O - electrons from Zakharov et al (2007)
!NUMBER OF COLL TRANS
%d
!NUMBER OF COLL TEMPS
%d
!COLL TEMPS""" % (freq.size, temp.size)
for i in temp: 
	print i,
print """
!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)"""
# deexcitation coefficients (i->j)
for i in range(freq.size):
	aij = pc.h*freq[i]/2./pc.k_B/temp
	gigj = weight[transitions[i][0]]/weight[transitions[i][1]]
	Ceij = ve*sigmaeij[i]*2.*aij*numpy.exp(aij)*special.k0(aij)
	Ceji = ve*gigj*sigmaeij[i]*2.*aij*numpy.exp(-aij)*special.k0(aij)
	print i+1, transitions[i][0]+1, transitions[i][1]+1,
        for j in Ceij: print "%.3e"%j,
        print
# 	print 2*i+2, transitions[i][1]+1, transitions[i][0]+1,
#         for j in Ceji: print "%.3e"%j,
#         print
