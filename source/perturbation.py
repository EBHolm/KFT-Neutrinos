import numpy as np 
from utils import sph2cart, get_angular_binning, get_y_binning
from model import *
import multiprocessing as mp
import itertools
import os
from model import Green

####################
#
#   First order computation
#
####################

def FermiDirac(y):
    return 1/(1 + np.exp(y))

def LapNFW(x, z):
    # Laplacian of the NFW potential
    x1, x2, x3 = x[0], x[1], x[2]
    r = np.sqrt(x1**2 + x2**2 + x3**2)
    return 4*np.pi*G*rho0(z)*Rs(z)**3/(r*(r + Rs(z))**2)/speedoflight**2*m_to_kpc

def get_first_order_correction(p_ini, x_here, z_span, N_z, z_binning='log', logz_min=-12):
    # Computes the first order correction to the density
    first_order = 0
    if z_binning == 'linear':
        zlist = np.linspace(z_span[0], z_span[1], N_z) # assume z-list is uniform 
        dz    = zlist[0] - zlist[1]
        first_order += 0.5*dz*mnu*(Green(0) - Green(zlist[0]))*LapNFW(x_here - (Green(0) - Green(zlist[0]))*p_ini, zlist[0])/(1 + zlist[0])/H(zlist[0])
        first_order += 0.5*dz*mnu*(Green(0) - Green(zlist[-1]))*LapNFW(x_here - (Green(0) - Green(zlist[-1]))*p_ini, zlist[-1])/(1 + zlist[-1])/H(zlist[-1])
        for z in zlist[1:-1]:
            first_order += dz*mnu*(Green(0) - Green(z))*LapNFW(x_here - (Green(0) - Green(z))*p_ini, z)/(1 + z)/H(z)
            
    elif z_binning == 'log':
        # assume lower bound on z_span is 0
        logzlist = np.linspace(np.log(np.max(z_span)), logz_min, N_z)
        zlist = np.exp(logzlist)
        dlogz = np.abs(logzlist[1] - logzlist[0])
        first_order += 0.5*dlogz*zlist[0]*mnu*(Green(0) - Green(zlist[0]))*LapNFW(x_here - (Green(0) - Green(zlist[0]))*p_ini, zlist[0])/(1 + zlist[0])/H(zlist[0])
        first_order += 0.5*dlogz*zlist[-1]*mnu*(Green(0) - Green(zlist[-1]))*LapNFW(x_here - (Green(0) - Green(zlist[-1]))*p_ini, zlist[-1])/(1 + zlist[-1])/H(zlist[-1])
        for z in zlist[1:-1]:
            first_order += dlogz*z*mnu*(Green(0) - Green(z))*LapNFW(x_here - (Green(0) - Green(z))*p_ini, z)/(1 + z)/H(z)
    return first_order

def get_first_order_integral(p_bins, T, x_here, z_span, N_z, y_binning='laggauss', y_max=10, z_binning='log', logz_min=-12):
    N_y, N_theta, N_phi = p_bins[0], p_bins[1], p_bins[2]
    y_list, y_weights = get_y_binning(N_y, y_max=y_max, method=y_binning)
    theta_list, theta_weights, phi_list, phi_weights = get_angular_binning(N_theta, N_phi)
    determinants = np.zeros([N_y, N_theta, N_phi])
    int_first = 0
    for idx_theta, (theta, theta_weight) in enumerate(zip(theta_list, theta_weights)):
        for idx_phi, (phi, phi_weight) in enumerate(zip(phi_list, phi_weights)):
            for idx_y, (y, y_weight) in enumerate(zip(y_list, y_weights)):
                p_ini = sph2cart(np.array([T*y, theta, phi]))
                weight = T**3*np.exp(y)*np.sin(theta)*y_weight*theta_weight*phi_weight
                first_order = get_first_order_correction(p_ini, x_here, z_span, N_z, z_binning, logz_min=logz_min)
                int_first += weight*y**2*FermiDirac(y)*(1 + first_order)
                determinants[idx_y, idx_theta, idx_phi] = 1/(1 + first_order)
    return int_first, determinants


####################
#
#   Second order computation
#
####################

def R1(r, z):
    front_factor = 4*np.pi*G*rho0(z)*Rs(z)**3/r**3/speedoflight**2*m_to_kpc
    return front_factor*(np.log((r + Rs(z))/Rs(z)) - r/(r + Rs(z)))

def R2(r, z):
    front_factor = 4*np.pi*G*rho0(z)*Rs(z)**3/r**3/speedoflight**2*m_to_kpc
    return front_factor*(3*(r + Rs(z))**2*np.log((r + Rs(z))/Rs(z)) - r*(4*r + 3*Rs(z)))/(r**2*(r + Rs(z))**2)

def R3(r, z):
    front_factor = 4*np.pi*G*rho0(z)*Rs(z)**3/r**3/speedoflight**2*m_to_kpc
    return front_factor*(15*(r + Rs(z))**3*np.log((r + Rs(z))/Rs(z)) - r*(23*r**2 + 36*r*Rs(z) + 15*Rs(z)**2))/(r**4*(r + Rs(z))**3)

def D1V(x, R1):
    return x*R1

def D2V(x, R1, R2):
    out = np.zeros([3, 3])
    for i in range(3):
        for j in range(3):
            out[i, j] = x[i]*x[j]*R2
            if i == j:
                out[i, j] += R1
    return out

def D3V(x, R2, R3):
    out = np.zeros([3, 3, 3])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                out[i, j, k] = x[i]*x[j]*x[k]*R3
                if i == j:
                    out[i, j, k] += x[k]*R2
                elif j == k:
                    out[i, j, k] += x[i]*R2
                elif k == i:
                    out[i, j, k] += x[j]*R2
    return out

def get_second_order_correction(p_ini, x_here, z_span, Nbins):
    zlist = np.linspace(z_span[0], z_span[1], Nbins)
    dz1 = dz2 = zlist[0] - zlist[1]
    second_order = 0
    for idx_2, z2 in enumerate(zlist):
        x2 = x_here - (Green(0) - Green(z2))*p_ini
        r2 = np.sqrt(np.dot(x2, x2))
        R2_x2 = R2(r2, z2)
        R3_x2 = R3(r2, z2)

        for idx_1, z1 in enumerate(zlist[0:idx_2]):
            x1 = x_here - (Green(0) - Green(z1))*p_ini
            r1 = np.sqrt(np.dot(x1, x1))
            R1_x1 = R1(r1, z1)
            
            first_term  = (Green(0) - Green(z1))*LapNFW(x2, z2)*LapNFW(x1, z1)
            second_term = (Green(z2) - Green(z1))*np.dot(np.trace(D3V(x2, R2_x2, R3_x2), axis1=1, axis2=2), D1V(x1, R1_x1))
            
            second_order += (first_term + second_term)*mnu/(1 + z1)/H(z1)*mnu/(1 + z2)/H(z2)*dz1*dz2
    return second_order

def second_order_integrand(inp):
    p_ini, x_here, z1, z2 = inp[0], inp[1], inp[2][0], inp[2][1]
    x2 = x_here - (Green(0) - Green(z2))*p_ini
    r2 = np.sqrt(np.dot(x2, x2))
    R2_x2 = R2(r2, z2)
    R3_x2 = R3(r2, z2)
    x1 = x_here - (Green(0) - Green(z1))*p_ini
    r1 = np.sqrt(np.dot(x1, x1))
    R1_x1 = R1(r1, z1)

    first_term  = (Green(0) - Green(z1))*LapNFW(x2, z2)*LapNFW(x1, z1)
    second_term = (Green(z2) - Green(z1))*np.dot(np.trace(D3V(x2, R2_x2, R3_x2), axis1=1, axis2=2), D1V(x1, R1_x1))

    integrand = (first_term + second_term)*mnu/(1 + z1)/H(z1)*mnu/(1 + z2)/H(z2)
    return integrand 

def get_second_order_correction_parallel(p_ini, x_here, z_span, Nbins, Ncores='all'):
    if Ncores == 'all':
        N_processes = os.cpu_count()
    elif Ncores == 'slurm':
        N_processes = os.environ['slurm_cpus_per_task']

    zlist = np.linspace(z_span[0], z_span[1], Nbins)
    dz1 = dz2 = zlist[0] - zlist[1]
    zpairs = []
    for idx_2, z2 in enumerate(zlist):
        for z1 in zlist[0:idx_2]:
            zpairs.append((z1, z2))
    inputs = [(p_ini, x_here, zs) for zs in zpairs]

    # if __name__ == '__main__': # Use this if running directly from this script
    if __name__ == 'relic_density':
        p = mp.Pool(N_processes)

        # Compute determinants, parallellized
        integrands = p.map(second_order_integrand, inputs)
        p.close()
        p.join()

        second_order = np.sum(np.array(integrands)*dz1*dz2)
        return second_order

def get_second_order_integral(p_bins, T, x_here, z_span, N_z, y_binning='laggauss', y_max=10):
    N_y, N_theta, N_phi = p_bins[0], p_bins[1], p_bins[2]
    y_list, y_weights = get_y_binning(N_y, y_max=y_max, method=y_binning)
    theta_list, theta_weights, phi_list, phi_weights = get_angular_binning(N_theta, N_phi)
    determinants = np.zeros([N_y, N_theta, N_phi])
    int_second = 0
    for idx_theta, (theta, theta_weight) in enumerate(zip(theta_list, theta_weights)):
        for idx_phi, (phi, phi_weight) in enumerate(zip(phi_list, phi_weights)):
            for idx_y, (y, y_weight) in enumerate(zip(y_list, y_weights)):
                p_ini = sph2cart(np.array([T*y, theta, phi]))
                weight = T**3*np.exp(y)*np.sin(theta)*y_weight*theta_weight*phi_weight
                first_order = get_first_order_correction(p_ini, x_here, z_span, N_z)
                second_order = get_second_order_correction_parallel(p_ini, x_here, z_span, N_z)
                int_second += weight*y**2*FermiDirac(y)*(1 + first_order + second_order)
                determinants[idx_y, idx_theta, idx_phi] = 1/(1 + first_order + second_order)
    return int_second, determinants

