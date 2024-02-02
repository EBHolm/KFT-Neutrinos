//
//  perturbation_second_kft.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 03/03/2023.
//

#ifndef perturbation_second_kft_hpp
#define perturbation_second_kft_hpp

#include <numbers>
#include "quadrature.hpp"

struct SecondOrderArgumentsKFT {
    double rtols[4];
    double atols[4];
    int GaussLaguerreNodes;
    
    double mnu;
    double Tnu;
    double r_here;
    double rr_here;
    double Mvir_over_Msun;
    int terms_flag;
    double G0;
    double z_ini;
    
    double weight;
    
    double Gz2;
    double G0_Gz2;
    double Rs2;
    double Rvir2;
    double front_factor2;
    
    double G0_Gz1;
    double Rs1;
    double Rvir1;
    double z1;
    double front_factor1; // = 4.*std::numbers::pi*rho0(z, Rs, conc)*pow(Rs, 3.)*pow(r, -3.);
    
    
    double p;
    double gp1;
    double gp2;
    double y_a1;
    double y_a2;
    
    std::vector<double> LaguerreNodes;
    std::vector<double> LaguerreWeights;
};


double second_order_kft(double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag = 0, double Tnu = 0.0001676375864435959);

double integrand_z2_kft(double z2, SecondOrderArgumentsKFT args);

double integrand_z1_kft(double z1, SecondOrderArgumentsKFT args);

double integrand_y_kft(double y, SecondOrderArgumentsKFT args);

double integrand_theta_kft(double theta, SecondOrderArgumentsKFT args);

double integrand_z2z1_kft(double z2, double z1, double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu = 0.0001676375864435959);

#endif /* perturbation_second_kft_hpp */
