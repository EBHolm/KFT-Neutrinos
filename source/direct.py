import multiprocessing as mp
import itertools
import os
import numpy as np 
from time import time 
from model import Green, derivatives
from utils import sph2cart, get_angular_binning, get_y_binning
import model 
Tnu = model.Tnu

from scipy.integrate import solve_ivp
def evolve_forward(x_ini, p_ini, z_span, abstol=1e-10, reltol=1e-12):
    # Solve Hamilton's equations
    initial = np.concatenate([x_ini, p_ini], axis=0)
    sol = solve_ivp(derivatives, z_span, initial, method='BDF', vectorized=True, atol=abstol, rtol=reltol)
    x_final = np.array([sol.y[i][-1] for i in [0, 1, 2]])
    p_final = np.array([sol.y[i][-1] for i in [3, 4, 5]])
    return x_final, p_final

def final_pos_residual(x_ini, p_ini, x_here, z_span, abstol=1e-12, reltol=1e-14):
    x_final, p_final = evolve_forward(x_ini, p_ini, z_span, abstol=abstol, reltol=reltol)
    return x_final - x_here
        
from scipy.optimize import root
def solve_boundary(p_ini, x_here, x_guess, z_span=[3, 0]):
    def root_func(x):
        return final_pos_residual(x, p_ini, x_here, z_span)
    return root(root_func, x_guess)

def get_determinant(x_ini, p_ini, eps, z_span, abstol=1e-10, reltol=1e-12):
    # Compute determinant in (4.8) with central differences
    derivs = np.zeros([3, 3])
    for i in range(3):
        # Vary initial position to compute derivative
        x_ini_var1, x_ini_var2 = np.copy(x_ini), np.copy(x_ini)
        x_ini_var1[i]   += eps
        x_ini_var2[i]   -= eps
        xf_var1, pf_var1 = evolve_forward(x_ini_var1, p_ini, z_span, abstol=abstol, reltol=reltol)
        xf_var2, pf_var2 = evolve_forward(x_ini_var2, p_ini, z_span, abstol=abstol, reltol=reltol)
        derivs[i,:]    = (xf_var1 - xf_var2)/(2*eps)
    return np.abs(np.linalg.det(derivs))

def det_from_p_ini(inputs):
    # Required for the parallellization in the next function
    p_ini = inputs[0]
    x_here, guess, eps, z_span, abstol, reltol, avoid_close = [inputs[i] for i in range(1, len(inputs))]
    p, theta, phi = p_ini[0], p_ini[1], p_ini[2]
    p_ini = sph2cart([p, theta, phi])

    # Compute a nice guess!
    guess_max = 1e+7
    if avoid_close: 
        # avoids guesses being close to halo. This usually happens for complex trajectories, and by excluding them, we help singling out the free branch
        free_guess = x_here - (Green(z_span[1]) - Green(z_span[0]))*sph2cart(p_ini)
        while True:
            if np.dot(guess, x_here) >= np.dot(free_guess, x_here):
                break
            else:
                guess = [(np.random.rand(1) - 0.5)*2*guess_max, (np.random.rand(1) - 0.5)*2*guess_max, (np.random.rand(1) - 0.5)*2*guess_max]
    #adjusted_free_guess = x_here - 2*(Green(z_span[1]) - Green(z_span[0]))*sph2cart(p_ini)

    x_ini = solve_boundary(p_ini, x_here, guess).x
    return get_determinant(x_ini, p_ini, eps, z_span, abstol=abstol, reltol=reltol)

def compute_integral_parallel(N_bins, T, x_here, guess, eps, z_span, Ncores='all', y_binning='laggauss', y_max=10, abstol=1e-10, reltol=1e-12, avoid_close=False):
    # Carries out integral (4.8) in parallel using Python's multiprocessing module
    # Should work on clusters out of the box!
    N_y, N_theta, N_phi = N_bins[0], N_bins[1], N_bins[2]

    y_list, y_weights = get_y_binning(N_y, y_max=y_max, method=y_binning)
    theta_list, theta_weights, phi_list, phi_weights = get_angular_binning(N_theta, N_phi)

    y_inis    = list(itertools.product(y_list, theta_list, phi_list))
    y_weights = list(itertools.product(y_weights, theta_weights, phi_weights))
    
    # We integrate over y, but the BVP solves using p
    p_inis    = list(itertools.product(y_list*T, theta_list, phi_list))
    det_inputs    = [(p_ini, x_here, guess, eps, z_span, abstol, reltol, avoid_close) for p_ini in p_inis]
    
    if Ncores == 'all':
        N_processes = os.cpu_count()
    elif Ncores == 'slurm':
        N_processes = os.environ['slurm_cpus_per_task']

    # if __name__ == '__main__': # Use this if running directly from this script
    if __name__ == 'relic_density':
        p = mp.Pool(N_processes)

        # Compute determinants, parallellized
        determinants = p.map(det_from_p_ini, det_inputs)
        p.close()
        p.join()

        # Compute integral; takes negligible time, so serial is OK
        int, int_free = 0, 0
        for idx, (y_ini, weight) in enumerate(zip(y_inis, y_weights)):
            y, theta, phi = y_ini[0], y_ini[1], y_ini[2]
            if y <= min(y_list):
                continue
            if y_binning == 'laggauss':
                total_weight = np.exp(y)*np.sin(theta)*weight[0]*weight[1]*weight[2]
            else:
                total_weight = np.sin(theta)*weight[0]*weight[1]*weight[2]
            int += y**2/(1 + np.exp(y))/determinants[idx]*total_weight
            int_free += y**2/(1 + np.exp(y))*total_weight
        int *= Tnu**3
        int_free *= Tnu**3
        return int, int_free, np.array(determinants).reshape(N_y, N_theta, N_phi)

from scipy.interpolate import CubicSpline
def compute_integral_serial(N_bins, T, x_here, eps, z_span, interp=False, y_binning='laggauss', y_max=10, abstol=1e-10, reltol=1e-12, silent=False):
    # Serial version of `compute_integral_parallel` which uses the prior BVP solution as guess for the next one 
    # if interp is True, will throw away dets above 1 and cubic spline the leftovers
    N_y, N_theta, N_phi = N_bins[0], N_bins[1], N_bins[2]
    y_list, y_weights = get_y_binning(N_y, y_max=y_max, method=y_binning)
    theta_list, theta_weights, phi_list, phi_weights = get_angular_binning(N_theta, N_phi)

    # Compute determinants
    determinants = np.zeros([N_y, N_theta, N_phi])
    initials = []
    for idx_theta, (theta, theta_weight) in enumerate(zip(theta_list, theta_weights)):
        initials_theta = []
        for idx_phi, (phi, phi_weight) in enumerate(zip(phi_list, phi_weights)):
            # Start from largest y-value with the free-streaming guess 
            y_list, y_weights = np.flip(y_list), np.flip(y_weights)
            guesses = []

            timer = time()
            p_ini_0 = sph2cart([y_list[0]*T, theta, phi])
            guesses.append(x_here - 2*(Green(z_span[1]) - Green(z_span[0]))*p_ini_0) # this should catch the free branch in most cases
            x_ini_0 = solve_boundary(p_ini_0, x_here, guesses[0]).x
            determinants[0, idx_theta, idx_phi] = get_determinant(x_ini_0, p_ini_0, eps, z_span, abstol=abstol, reltol=reltol)
            if not silent:
                print(f"Iteration idx_y=0/{N_y-1}, idx_theta={idx_theta}/{N_theta-1}, idx_phi={idx_phi}/{N_phi-1}, det={determinants[0, idx_theta, idx_phi]:.3} (time elapsed: {time() - timer:.3} s).")
            for idx_y, y in enumerate(y_list[1:]):
                p_ini = sph2cart([y*T, theta, phi])
                x_ini = solve_boundary(p_ini, x_here, guesses[idx_y]).x
                guesses.append(x_ini)
                determinants[idx_y + 1, idx_theta, idx_phi] = get_determinant(x_ini, p_ini, eps, z_span, abstol=abstol, reltol=reltol)
                if not silent:
                    print(f"Iteration idx_y={idx_y+1}/{N_y-1}, idx_theta={idx_theta}/{N_theta-1}, idx_phi={idx_phi}/{N_phi-1}, det={determinants[idx_y + 1, idx_theta, idx_phi]:.3} (time elapsed: {time() - timer:.3} s).")
            initials_theta.append(guesses) # store initial guesses

            if interp:
                # This only really works if outliers appear in "holes" and have good points to either side
                
                # If the smallest y-point is divergent, set it equal to free value, since we anyway only have a lower bound
                if np.abs(determinants[-1, idx_theta, idx_phi]) > 1.0:
                    determinants[-1, idx_theta, idx_phi] = 1.0 # required for interpolating

                # Define outliers as being above 1.0, the free solution
                dets = np.flip(determinants[:, idx_theta, idx_phi])
                not_outlier_indices = dets < 1.0

                # Spline the non-outlier determinants
                y_pos = np.flip(y_list)
                cs = CubicSpline(y_pos[not_outlier_indices], dets[not_outlier_indices], bc_type='clamped')
                determinants[:, idx_theta, idx_phi] = cs(y_list)
        initials.append(initials_theta)

    # Carry out integration
    int, int_free = 0, 0
    for idx_theta, (theta, theta_weight) in enumerate(zip(theta_list, theta_weights)):
        for idx_phi, (phi, phi_weight) in enumerate(zip(phi_list, phi_weights)):
            for idx_y, (y, y_weight) in enumerate(zip(y_list, y_weights)):
                if y_binning == 'laggauss':
                    total_weight = np.exp(y)*np.sin(theta)*y_weight*theta_weight*phi_weight
                else:
                    total_weight = np.sin(theta)*y_weight*theta_weight*phi_weight
                int += y**2/(1 + np.exp(y))/determinants[idx_y, idx_theta, idx_phi]*total_weight
                int_free += y**2/(1 + np.exp(y))*total_weight
    int *= Tnu**3
    int_free *= Tnu**3
    return int, int_free, determinants[::-1, :, :], np.array(initials)
