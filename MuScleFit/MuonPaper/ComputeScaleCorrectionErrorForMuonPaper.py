# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 17:08:07 2011

@author: M. De Mattia - marco.de.mattia@cern.ch
"""

import math

b = -6.25*(10**(-5))
sigma_b = 0.24*(10**(-5))
c = 1.81*(10**(-3))
sigma_c = 0.21*(10**(-3))
d = 1.277*(10**(-4))
sigma_d = 0.050*(10**(-4))
e = 2.71*(10**(-1))
sigma_e = 0.25*(10**(-1))

pt = 45.6
pt = 30
eta = 2.1

# print "pi =", math.pi

sigma_f = pt*math.sqrt(sigma_b**2 + math.pi*sigma_d**2 +
                       math.pi*(d**2)*sigma_e**2 +
                       (eta**4)*(sigma_c**2)/(pt**2))

print "sigma_b**2 =", sigma_b**2
print "math.pi*sigma_d**2 =", math.pi*sigma_d**2
print "math.pi*(d**2)*sigma_e**2 =", math.pi*(d**2)*sigma_e**2
print "(eta**4)*(sigma_c**2)/(pt**2) =", (eta**4)*(sigma_c**2)/(pt**2)

print ""
print sigma_d*(pt**2)
print ""

print "deltaPt_amplitude =", d*pt**2
print "deltaPt(no phi) =", b*(pt**2) + c*pt*(eta**2)
# print "sigma_f_overPt =", sigma_f
print "sigma_f =", sigma_f*pt