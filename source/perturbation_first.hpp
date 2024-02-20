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
#include "cosmology.hpp"


struct FirstOrderArguments { // Could make this a struct now
    double rtols[3];
    double atols[3];
    int GaussLaguerreNodes;
    
    double mnu;
    double Tnu;
    double r_here;
    double rr_here;
    double Mvir_over_Msun;
    double conc;
    double Rvir;
    double front_factor;
    double G0_Gz;
    double weight;
    double Rs;
    double p;
    double gp;
    double r_a;
    double r;
    double phi_integrand;
    double theta;
    double z_ini;
};

double zeroth_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_y_zeroth(double y, FirstOrderArguments args);
double integrand_theta_zeroth(double theta, FirstOrderArguments args);

double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_z(double z, FirstOrderArguments args);
double integrand_y(double y, FirstOrderArguments args);
double integrand_theta(double theta, FirstOrderArguments args);
double integrand_phi(double phi, FirstOrderArguments args);

double epsilon_first(double y, double theta, FirstOrderArguments args);
double epsilon_first(double y, double theta, double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_z_epsilon(double z, FirstOrderArguments args);

double integrand_y_complete(double y, double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);
double integrand_z_complete(double z, FirstOrderArguments args);

double integrand_stub(double y, FirstOrderArguments args);

#endif /* perturbation_first_hpp */
