import numpy as np
H0 = 67.36; # Hubble constant in km/s/Mpc or m/s/kpc 
OmegaL = 0.6847
OmegaM = 0.3153
#Mvir = 3.76e12*1.988e30 # in units of kg 
Msun = 1.98847e+30 # Solar mass in kg
Mvir = 2.03*10**12*Msun
Mvir_over_Msun = 2.03*10**12

Tnu = 2.72548*pow(4./11.,1./3.)*8.61732814974056e-5 # relic neutrino temperature today in eV 
speedoflight = 2.99792458e8 # speed of light in m/s 
H0_kpc = H0/speedoflight
little_h = H0/100
G = 6.674e-11 # Gravitational constant in m^3/kg/s^2 
eV_to_inv_m = 5.06e+5
m_to_kpc = 3.24e-20

eV_to_inv_kpc = 2.489e+25

#mnu = 0.3 # neutrino mass in eV
mnu = 0.15
ToverM = Tnu/mnu
#mnu *= eV_to_inv_m/m_to_kpc # now in 1/kpc
#Tnu *= eV_to_inv_m/m_to_kpc


#####################################################
#
#   Problem definition functions
#   Captures the redshift dependence of the gravitational environment
#
#####################################################

def H(z):
    # Hubble expansion rate in [1/kpc]
    H = H0*np.sqrt(OmegaM*(1+z)**3+OmegaL)/speedoflight
    return H

def conc_old(z):
    # Concentration parameter in the NFW profile, (4.5) in Ringwald & Wong 2004
    c_0 = 9*(Mvir/(1.5e+13/(H0/100)*Msun))**(-0.13)
    return c_0/(1 + z)

def a_factor(z):
    a = 0.537 + (1.025 - 0.537)*np.exp(-0.718*z**(1.08))
    return 0.537 + (1.025 - 0.537)*np.exp(-0.718*z**1.08)

def b_factor(z):
    return -0.97 + 0.024*z

beta = 0.613
def conc(z):
    return np.exp(a_factor(z) + b_factor(z)*np.log(Mvir_over_Msun/(10**12*little_h**(-1))))/beta

def Omega_m(z):
    return OmegaM*(1 + z)**3/(OmegaM*(1 + z)**3 + OmegaL)

def rho_crit(z): # in kg/kpc^3
    return 3*H(z)**2/(8*np.pi*G*speedoflight**(-2)*m_to_kpc)

def Delta_vir(z): # unitless
    return 18*np.pi**2 + 82*(Omega_m(z) - 1) - 39*(Omega_m(z) - 1)**2

def Rvir(z):
    # virial radius in kpc
    return (3*Mvir/(4*np.pi*Delta_vir(z)*rho_crit(z)))**(1/3)

def Rs(z):
    # characteristic inner radius in kpc
    return Rvir(z)/conc(z)

def rho0(z): # in kg/kpc^3
    return Mvir*(4*np.pi*Rs(z)**3/(1 + z)**3*(np.log(1 + conc(z)) - conc(z)/(1 + conc(z))))**(-1)

def NFW(x, z):
    # Gradient of NFW potential; input always Cartesian
    x = np.array(x)
    r = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    return -x*4*np.pi*G*rho0(z)*Rs(z)**3/r**3*(r/(r + Rs(z)) - np.log(1 + r/Rs(z)))/speedoflight**2*m_to_kpc

def derivatives_free(z, y):
    x = np.array([y[0], y[1], y[2]])
    p = np.array([y[3], y[4], y[5]])
    dxdz = -(1 + z)*p/H(z)/mnu
    dpdz = np.zeros(3)
    return [*dxdz, *dpdz]
    
def derivatives(z, y):
    a = 1/(1 + z)
    x = np.array([y[0], y[1], y[2]])
    p = np.array([y[3], y[4], y[5]])
    dxdz = -(1 + z)*p/H(z)/mnu
    dpdz = mnu*NFW(x, z)/H(z)/(1 + z)
    return [*dxdz, *dpdz]

from scipy.special import hyp2f1
def Green(z):
    # Green's function for the free particle
    return -(1 + z)**2/(2*H0_kpc*mnu*np.sqrt(OmegaL))*hyp2f1(1/2, 2/3, 5/3, -OmegaM/OmegaL*(1 + z)**3)