import numpy as np 
from model import *
from quadrature import GaussKronrod_adapt, GaussLaguerre
from utils import sph2cart_quick
from collections import namedtuple

####################
#
#   First order computation
#
####################

# Integral is computed in the following sequence: z -> p -> theta -> phi
# Currently very hardcoded to NFW. How to generalize nicely?

KftArgs = namedtuple('KftArgs', 
                     field_names=['rtols', 'atols', 'f_evals', 'x_here', 'mnu', 'Tnu', 'G0_Gz', 
                                  'weight', 'Rs', 'p', 'theta', 'sintheta', 'costheta', 'zlist', 'Gauss_Laguerre'])

def integral(z_span, rtols, atols, x_here, mnu, Tnu, Gauss_Laguerre=0):
    a, b = z_span[0], z_span[1]
    args = KftArgs(rtols=rtols, atols=atols, x_here=x_here, mnu=mnu, Tnu=Tnu, f_evals=[0], zlist=[],
                G0_Gz=1, weight=1, Rs=1, p=1, theta=1, sintheta=1, costheta=1, Gauss_Laguerre=Gauss_Laguerre) # dummies
    assert(args)
    if Gauss_Laguerre:
        #I, err = GaussKronrod_adapt(integrand_z_GaussLaguerre, a, b, args.rtols[0], args.atols[0], is_indefinite=False, f_args=args)
        I, err = GaussKronrod_adapt(integrand_z_GaussLaguerre, a, b, args.rtols[0], args.atols[0], is_indefinite=False, f_args=args)
    else:
        I, err = GaussKronrod_adapt(integrand_z, a, b, args.rtols[0], args.atols[0], is_indefinite=False, f_args=args)
    return I, args.f_evals[0], args.zlist

def integrand_z(z, f_args):
    f_args.zlist.append(z)
    Hz = H(z)
    con = conc(z)
    Rs = Rs_quick(z, con, Hz)
    G0_Gz = Green_quick(0, f_args.mnu) - Green_quick(z, f_args.mnu)
    z_weight = f_args.mnu*G0_Gz/(1 + z)/Hz*4*np.pi*G*rho0_quick(z, Rs, con)*Rs**3/speedoflight**2*m_to_kpc
    
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=G0_Gz, weight=z_weight, Rs=Rs, zlist=f_args.zlist, Gauss_Laguerre=f_args.Gauss_Laguerre,
                   p=1, theta=1, sintheta=1, costheta=1) # dummies
    pmax = f_args.Tnu*np.log(1/(1e-10) - 1)*0.5 # factor 0.5 for precision
    I, err = GaussKronrod_adapt(integrand_p, *[0, pmax], f_args.rtols[1], f_args.atols[1], is_indefinite=False, f_args=args)
    #I, err = GaussKronrod_adapt(integrand_p, [0, 1], f_args.rtols[1], f_args.atols[1], is_indefinite=True, f_args=args)
    # print(f"{z=}, {I=}", flush=True)
    return I

def integrand_z_GaussLaguerre(z, f_args):
    f_args.zlist.append(z)
    Hz = H(z)
    con = conc(z)
    Rs = Rs_quick(z, con, Hz)
    G0_Gz = Green_quick(0, f_args.mnu) - Green_quick(z, f_args.mnu)
    """print(f"{z=}")
    print(f"{G0_Gz=}")
    print(f"{f_args.mnu=}")
    print(f"{Green(0.2, 0.1)=}")
    asdfasdf """
    z_weight = f_args.mnu*G0_Gz/(1 + z)/Hz*4*np.pi*G*rho0_quick(z, Rs, con)*Rs**3/speedoflight**2*m_to_kpc
    
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=G0_Gz, weight=z_weight, Rs=Rs, zlist=f_args.zlist, Gauss_Laguerre=f_args.Gauss_Laguerre,
                   p=1, theta=1, sintheta=1, costheta=1) # dummies
    
    assert(type(f_args.Gauss_Laguerre) == int)
    assert(f_args.Gauss_Laguerre > 0)
    I = GaussLaguerre(integrand_y, N=f_args.Gauss_Laguerre, f_args=args)
    #print(f"{z=}, {I=}", flush=True)
    #asdf
    return I

def integrand_p(p, f_args):
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=f_args.G0_Gz, weight=f_args.weight*p**2/(1 + np.exp(p/f_args.Tnu)), Rs=f_args.Rs, p=p, zlist=f_args.zlist,
                   theta=1, sintheta=1, costheta=1, Gauss_Laguerre=f_args.Gauss_Laguerre) # dummies
    I, err = GaussKronrod_adapt(integrand_theta, 0, np.pi, f_args.rtols[2], f_args.atols[2], is_indefinite=False, f_args=args)
    # print(f"{p=}, {I=}", flush=True)
    return I

def integrand_y(y, f_args):
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=f_args.G0_Gz, weight=f_args.Tnu**3*f_args.weight*np.exp(y)*y**2/(1 + np.exp(y)), Rs=f_args.Rs, p=f_args.Tnu*y, zlist=f_args.zlist,
                   theta=1, sintheta=1, costheta=1, Gauss_Laguerre=f_args.Gauss_Laguerre) # dummies
    I, err = GaussKronrod_adapt(integrand_theta, 0, np.pi, f_args.rtols[2], f_args.atols[2], is_indefinite=False, f_args=args)
    # print(f"{y=}, {I=}", flush=True)
    return I

def integrand_theta(theta, f_args):
    sintheta = np.sin(theta)
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=f_args.G0_Gz, weight=f_args.weight*sintheta, Rs=f_args.Rs, p=f_args.p, zlist=f_args.zlist, Gauss_Laguerre=f_args.Gauss_Laguerre,
                   theta=theta, sintheta=sintheta, costheta=np.cos(theta))
    #print(f"weight = {f_args.weight}")
    #asdf
    I, err = GaussKronrod_adapt(integrand_phi, 0, 2*np.pi, f_args.rtols[3], f_args.atols[3], is_indefinite=False, f_args=args)
    # print(f"{theta=}, {I=}")
    return I

def integrand_phi(phi, f_args):
    f_args.f_evals[0] += 1
    p_ini = sph2cart_quick([f_args.p, f_args.theta, phi], f_args.sintheta, f_args.costheta)
    pos_back = f_args.x_here - f_args.G0_Gz*p_ini
    print(f"{p_ini=}")
    print(f"{phi=}")
    print(f"{f_args.G0_Gz=}")
    print(f"pos_back = {pos_back}")
    print(f"r = {np.sqrt(pos_back[0]**2 + pos_back[1]**2 + pos_back[2]**2)}")
    print(f"Rs = {f_args.Rs}")
    print(f"weight = {f_args.weight}")
    print(f"{LapNFW_quick(f_args.x_here - f_args.G0_Gz*p_ini, f_args.Rs)=}")
    asdf
    return f_args.weight*LapNFW_quick(f_args.x_here - f_args.G0_Gz*p_ini, f_args.Rs)

def integrand_stub(x, f_args):
    return 1.0

#########################################################

"""

KftArgs = namedtuple('KftArgs', 
                     field_names=['rtols', 'atols', 'f_evals', 'x_here', 'mnu', 'Tnu', 'G0_Gz', 
                                  'weight', 'Rs', 'p', 'theta', 'sintheta', 'costheta'])

def integral(z_span, rtols, atols, x_here, mnu, Tnu):
    args = KftArgs(rtols=rtols, atols=atols, x_here=x_here, mnu=mnu, Tnu=Tnu, f_evals=[0],
                   G0_Gz=1, weight=1, Rs=1, p=1, theta=1, sintheta=1, costheta=1) # dummies
    I, err = GaussKronrod_adapt(integrand_z, z_span, args.rtols[0], args.atols[0], is_indefinite=False, f_args=args)
    return 1+I, args.f_evals[0]

def integrand_z(z, f_args):
    Hz = H(z)
    con = conc(z)
    Rs = Rs_quick(z, con, Hz)
    G0_Gz = Green_quick(0, f_args.mnu) - Green_quick(z, f_args.mnu)
    z_weight = f_args.mnu*G0_Gz/(1 + z)/Hz*4*np.pi*G*rho0_quick(z, Rs, con)*Rs**3/speedoflight**2*m_to_kpc
    
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=G0_Gz, weight=z_weight, Rs=Rs,
                   p=1, theta=1, sintheta=1, costheta=1) # dummies
    I, err = GaussKronrod_adapt(integrand_p, [0, 1], f_args.rtols[1], f_args.atols[1], is_indefinite=True, f_args=args)
    print(f"{z=}, {I=}", flush=True)
    return I

def integrand_p(p, f_args):
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=f_args.G0_Gz, weight=f_args.weight*p**2/(1 + np.exp(p/f_args.Tnu)), Rs=f_args.Rs, p=p, 
                   theta=1, sintheta=1, costheta=1) # dummies
    I, err = GaussKronrod_adapt(integrand_theta, [0, np.pi], f_args.rtols[2], f_args.atols[2], is_indefinite=False, f_args=args)
   # print(f"{p=}, {I=}", flush=True)
    return I

def integrand_theta(theta, f_args):
    sintheta = np.sin(theta)
    args = KftArgs(rtols=f_args.rtols, atols=f_args.atols, x_here=f_args.x_here, mnu=f_args.mnu, Tnu=f_args.Tnu, f_evals=f_args.f_evals,
                   G0_Gz=f_args.G0_Gz, weight=f_args.weight*sintheta, Rs=f_args.Rs, p=f_args.p, 
                   theta=theta, sintheta=sintheta, costheta=np.cos(theta))
    I, err = GaussKronrod_adapt(integrand_phi, [0, 2*np.pi], f_args.rtols[3], f_args.atols[3], is_indefinite=False, f_args=args)
    #print(f"{theta=}, {I=}")
    return I

def integrand_phi(phi, f_args):
    f_args.f_evals[0] += 1
    p_ini = sph2cart_quick([f_args.p, f_args.theta, phi], f_args.sintheta, f_args.costheta)
    return f_args.weight*LapNFW_quick(f_args.x_here - f_args.G0_Gz*p_ini, f_args.Rs)

"""