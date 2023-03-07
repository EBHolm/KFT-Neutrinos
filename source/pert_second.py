import numpy as np 
from model import *
from quadrature import GaussKronrod_adapt
from utils import sph2cart

####################
#
#   Second order computation
#
####################

# Inspired by the first order computation

def integral(z_span, fargs):
    f_evals = 0
    I, err = GaussKronrod_adapt(integrand_z2, z_span, fargs[0][0], is_indefinite=False, f_args=fargs)
    return I, f_evals

def integrand_z2(z2, fargs):
    # tols = fargs[0]
    # x_here = fargs[1]
    # mnu = fargs[2]
    # Tnu = fargs[3]

    G0 = Green(0)
    G0_Gz2 = G0 - Green(z2)

    # z2-dependent potential quantities...

    fargs.append(G0_Gz2*fargs[2]/(1 + z2)/H(z2)) # z2_weight
    fargs.append(G0)
    fargs.append(G0_Gz2)
    I, err = GaussKronrod_adapt(integrand_z1, [0, 1], fargs[0][1], is_indefinite=False, f_args=fargs)
    return I

def integrand_z1(z1, fargs):
    fargs[4] *= fargs[2]/(1 + z1)/H(z1) # z1_weight

    # z1-dependent potential quantities...

    fargs.append(fargs[5] - Green(z1)) # G0_Gz1
    I, err = GaussKronrod_adapt(integrand_p, [0, 1], fargs[0][2], is_indefinite=True, f_args=fargs)
    return I

def integrand_p(p, fargs):
    fargs[4] *= p**2/(1 + np.exp(p/fargs[3])) # p_weight(p)
    fargs.append(p)
    I, err = GaussKronrod_adapt(integrand_theta, [0, np.pi], fargs[0][3], is_indefinite=False, f_args=fargs)
    return I

def integrand_theta(theta, fargs):
    fargs[4] *= np.sin(theta) # theta_weight(theta)
    fargs.append(theta)
    I, err = GaussKronrod_adapt(integrand_phi, [0, 2*np.pi], fargs[0][4], is_indefinite=False, f_args=fargs)
    return I

def integrand_phi(phi, fargs):
    # weight = fargs[4]
    # G0 = fargs[5]
    # G0-Gz2 = fargs[6]
    # G0-Gz1 = fargs[7]
    # p = fargs[8]
    # theta = fargs[9]
    p_ini = sph2cart(np.array([fargs[7], fargs[8], phi]))

    return 0.0


