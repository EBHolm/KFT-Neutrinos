import numpy as np 

# Coordinate transformations
def sph2cart(x):
    r, th, phi = x[0], x[1], x[2]
    x1 = r*np.sin(th)*np.cos(phi)
    x2 = r*np.sin(th)*np.sin(phi)
    x3 = r*np.cos(th)
    return np.array([x1, x2, x3])

def cart2sph(x):
    x1, x2, x3 = x[0], x[1], x[2]
    r   = np.sqrt(x1**2 + x2**2 + x3**2)
    th  = np.arctan(np.sqrt(x1**2 + x2**2)/x3)
    phi = np.arctan(x2/x3)
    return np.array([r, th, phi])


from numpy.polynomial.laguerre import laggauss
def get_y_binning(N_y, y_max=10, method='laggauss'):
    if method == 'laggauss':
        y_list, y_weights = laggauss(N_y)
    elif method == 'linear':
        #y_list        = np.array([j*y_max/N_y for j in range(0, N_y)])
        #y_weights     = y_max/N_y*np.ones(N_y)
        y_list     = np.linspace(y_max/N_y, y_max, N_y)
        y_weights  = np.ones(N_y)*(y_list[1] - y_list[0])
        y_weights[0]  /= 2
        y_weights[-1] /= 2
    elif method == 'log':
        y_list = np.logspace(-1.5, np.log(y_max), N_y, endpoint=False)
        y_weights = np.diff(y_list)
        y_weights = np.append(y_weights, 0) # hopefully, y_max is large enough that the last bin doesn't contribute anyway 
    return y_list, y_weights

def get_angular_binning(N_theta, N_phi):
    # k=0 and k=Ntheta give sin(theta)=0, so maybe omit these? Well, the determinant can make up for this cancellation in principle
    # TODO: Implement Simpson's rule. Could also consider making sin(theta) into a quadrature weight!
    if N_theta > 1 and N_phi > 1:
        theta_list    = np.linspace(0, np.pi, N_theta, endpoint=True)
        phi_list      = np.linspace(0, 2*np.pi, N_phi, endpoint=True)
        theta_weights = np.ones(N_theta)*(theta_list[1] - theta_list[0]) # uniform grids 
        phi_weights   = np.ones(N_phi)*(phi_list[1] - phi_list[0])
        theta_weights[0]  /= 2 # trapezoidal rule
        phi_weights[0]    /= 2
        theta_weights[-1] /= 2
        phi_weights[-1]   /= 2
    else: # only one bin requested
        if N_theta == 1 and N_phi > 1:
            theta_list = [np.pi/2]
            theta_weights = [np.pi]
            phi_list      = np.linspace(0, 2*np.pi, N_phi, endpoint=True)
            phi_weights   = np.ones(N_phi)*(phi_list[1] - phi_list[0])
            phi_weights[0]  /= 2
            phi_weights[-1] /= 2
        elif N_phi == 1 and N_theta > 1:
            phi_weights = [2*np.pi]
            phi_list = [np.pi]
            theta_list    = np.linspace(0, np.pi, N_theta, endpoint=True)
            theta_weights = np.ones(N_theta)*(theta_list[1] - theta_list[0]) # uniform grids 
            theta_weights[0]  /= 2
            theta_weights[-1] /= 2
        elif N_phi == 1 and N_theta == 1:
            theta_list = [np.pi/2]
            #theta_list = [0.0]
            phi_list = [np.pi]
            #phi_list = [0.0]
            theta_weights = [np.pi]
            phi_weights = [2*np.pi]
    return theta_list, theta_weights, phi_list, phi_weights
