import cython
import numpy as np 

cdef extern from "perturbation_first.hpp":
    cdef double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, int N_GaussLaguerre, double Tnu)
    
cpdef py_first_order(double mass, double z_ini, double[::1] rtols, double[::1] atols, double r_here, int N_GaussLaguerre, double Tnu=0.0001676375864435959): 
    integral = first_order(mass, z_ini, &rtols[0], &atols[0], r_here, N_GaussLaguerre, Tnu)
    return integral
