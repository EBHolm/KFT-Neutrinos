//
//  perturbation_first.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#ifndef perturbation_first_hpp
#define perturbation_first_hpp

#include <numbers>

#include "quadrature.hpp"
#include "utils.hpp"
#include "cosmology.hpp"


struct FirstOrderArguments { // Could make this a struct now
    double rtols[3];
    double atols[3];
    int GaussLaguerreNodes;
    
    double mnu;
    double Tnu;
    double r_here;
    double rr_here;
    
    double conc;
    double Rvir;
    double z_factor;
    double G0_Gz;
    double weight;
    double Rs;
    double p;
    double gp;
    double r_a;
    double r;
    double phi_integrand;
};

double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_z(double z, FirstOrderArguments args);
double integrand_y(double y, FirstOrderArguments args);
double integrand_theta(double theta, FirstOrderArguments args);
double integrand_phi(double phi, FirstOrderArguments args);

double integrand_y_complete(double y, double mass, double z_ini, double rtols[3], double atols[3], double r_here, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_z_complete(double z, FirstOrderArguments args);

double integrand_stub(double y, FirstOrderArguments args);

#endif /* perturbation_first_hpp */
