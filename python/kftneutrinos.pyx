import cython
import numpy as np 

cdef extern from "perturbation_first.hpp":
    cdef double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu)
    cdef double integrand_y_complete(double y, double mass, double z_ini, double rtols[2], double atols[2], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu);

cdef extern from "perturbation_second_kft.hpp":
    cdef double second_order_kft(double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag, double Tnu)
    cdef double integrand_z2z1_kft(double z2, double z1, double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu)

cpdef py_first_order(double mass, double z_ini, double[::1] rtols, double[::1] atols, double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu=0.0001676375864435959): 
    integral = first_order(mass, z_ini, &rtols[0], &atols[0], r_here, Mvir_over_Msun, N_GaussLaguerre, Tnu)
    return integral

cpdef py_first_integrand_y(double y, double mass, double z_ini, double[::1] rtols, double[::1] atols, double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu=0.0001676375864435959):
    integrand = integrand_y_complete(y, mass, z_ini, &rtols[0], &atols[0], r_here, Mvir_over_Msun, N_GaussLaguerre, Tnu)
    return integrand

cpdef py_second_order_kft(double mass, double z_ini, double[::1] rtols, double[::1] atols, double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag=0, double Tnu=0.0001676375864435959):
    integral = second_order_kft(mass, z_ini, &rtols[0], &atols[0], r_here, Mvir_over_Msun, N_GaussLaguerre, terms_flag, Tnu)
    return integral

cpdef py_second_integrand_z1z2_kft(double z2, double z1, double mass, double z_ini, double[::1] rtols, double[::1] atols, double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu=0.0001676375864435959):
    integrand = integrand_z2z1_kft(z2, z1, mass, z_ini, &rtols[0], &atols[0], r_here, Mvir_over_Msun, N_GaussLaguerre, Tnu)
    return integrand
