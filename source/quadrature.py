import numpy as np 
from numpy.polynomial.laguerre import laggauss

def GaussKronrod_quad(f, a, b, is_indefinite, f_args):
    z_k = np.array([-0.991455371120813, -0.949107912342759, -0.864864423359769, -0.741531185599394, 
                    -0.586087235467691, -0.405845151377397, -0.207784955007898, 0.0,
                    0.207784955007898, 0.405845151377397, 0.586087235467691, 0.741531185599394,
                    0.864864423359769, 0.949107912342759, 0.991455371120813])
    w_k = np.array([0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525,
                    0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728,
                    0.204432940075298, 0.190350578064785, 0.169004726639267, 0.140653259715525,
                    0.104790010322250, 0.063092092629979, 0.022935322010529])
    w_g = np.array([0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469,
                    0.381830050505119, 0.279705391489277, 0.129484966168870])
    
    # Transform z into t in interval between a and b:
    t = 0.5*(a*(1 - z_k) + b*(1 + z_k))
    # Modify weights such that it reflects the linear transformation above: */
    w_k *= 0.5*(b - a)
    w_g *= 0.5*(b - a)
    if is_indefinite:
        # Transform t into x in interval between 0 and inf:
        x = 1.0/t - 1.0
        # Modify weight accordingly:
        w_k /= (t*t)
        w_g /= (t[1::2]*t[1::2])
    else:
        x = t
    #assert(f_args)
    y = np.array([f(xx, f_args=f_args) for xx in x])
    
    Ik = np.dot(y, w_k)
    Ig = np.dot(y[1::2], w_g)
    err = (200*abs(Ik - Ig))**1.5
    return Ik, err

def GaussKronrod_adapt(f, a, b, rtol, abstol, f_args, is_indefinite=False):
    # Do adaptive Gauss-Kronrod quadrature.
    # To integrate between 0 and infinity, set a=0, b=1 and is_indefinite=True on first call
    I, err = GaussKronrod_quad(f, a, b, is_indefinite, f_args)
    # Stop recursion if tolerance has been made or recursion is too deep:
    
    if (abs(I) > 0 and err/I < rtol) or err < abstol:
        return I, err
    else:
        m = 0.5*(a + b)
        I_left, err_left = GaussKronrod_adapt(f, a, m, 1.5*rtol, is_indefinite, f_args)
        I_rght, err_rght = GaussKronrod_adapt(f, m, b, 1.5*rtol, is_indefinite, f_args)
        return I_left + I_rght, np.sqrt(err_left**2 + err_rght**2)
    
def GaussLaguerre(f, N, f_args):
    x_list, x_weights = laggauss(N)
    val = 0.0
    for x, weight in zip(x_list, x_weights):
        val += f(x, f_args)*weight
    return val