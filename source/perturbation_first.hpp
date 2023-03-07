//
//  perturbation_first.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#ifndef perturbation_first_hpp
#define perturbation_first_hpp

#include <cmath>
#include <functional>
#include <algorithm>
#include <numbers>

#include "quadrature.hpp"
#include "utils.hpp"
#include "cosmology.hpp"


double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_z(double z, FirstOrderArguments args);
double integrand_y(double y, FirstOrderArguments args);
double integrand_theta(double theta, FirstOrderArguments args);
double integrand_phi(double phi, FirstOrderArguments args);

double integrand_stub(double y, FirstOrderArguments args);

#endif /* perturbation_first_hpp */
