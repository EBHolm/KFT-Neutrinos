//
//  perturbation_second.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 03/03/2023.
//

#ifndef perturbation_second_hpp
#define perturbation_second_hpp

#include <stdio.h>
#include <numbers>
#include "quadrature.hpp"


double second_order(double mass, double z_ini, double rtols[4], double atols[4], double r_here, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);

double integrand_z2(double z2, SecondOrderArguments args);

double integrand_z1(double z1, SecondOrderArguments args);

double integrand_y(double y, SecondOrderArguments args);

double integrand_theta(double theta, SecondOrderArguments args);

double integrand_phi(double phi, SecondOrderArguments args);



#endif /* perturbation_second_hpp */
