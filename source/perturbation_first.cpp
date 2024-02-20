//
//  perturbation_first.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include "perturbation_first.hpp"

double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu) {
    FirstOrderArguments args = {{rtols[0], rtols[1], rtols[2]},
                                {atols[0], atols[1], atols[2]},
                                N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun};
    
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_z, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
};

double integrand_z(double z, FirstOrderArguments args) {
    double Hz = H(z);
    double con = conc(z, args.Mvir_over_Msun);
    args.conc = con;
    args.G0_Gz = Green(0, args.mnu) - Green(z, args.mnu);
    args.Rs = Rs(z, con, Hz, args.Mvir_over_Msun);
    args.Rvir = args.Rs*args.conc;
    args.front_factor = 4.0*std::numbers::pi*_G_*rho0(z, args.Rs, con, args.Mvir_over_Msun)*pow(args.Rs, 3.)*pow(1 + z, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight = args.mnu*args.G0_Gz/Hz/(1 + z);
    // double I = GaussLaguerre<FirstOrderArguments>(integrand_y, args.GaussLaguerreNodes, args);
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_y, 1e-8, 20.0, args.rtols[0], args.atols[0], args);
    return I;
};

double integrand_y(double y, FirstOrderArguments args) {
    args.p = args.Tnu*y;
    args.gp = args.p*args.G0_Gz;
    args.r_a = args.rr_here + pow(args.gp, 2.);
    args.weight *= pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    // args = first_update_args_y(y, args);
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_theta, 0.0, std::numbers::pi, args.rtols[1], args.atols[1], args);
    return I;
};

double integrand_theta(double theta, FirstOrderArguments args) {
    args.r = sqrt(args.r_a - 2.*args.gp*args.r_here*cos(theta));
    args.weight *= sin(theta);
    // return 2.*std::numbers::pi*args.weight*args.front_factor*LapNFWKepler(args.r, args.Rs, args.Rvir);
    return 2.*std::numbers::pi*args.weight*args.front_factor*LapNFW(args.r, args.Rs);
};

double integrand_y_complete(double y, double mass, double z_ini, double rtols[2], double atols[2], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu) {
    /* Returns d delta n/dp, i.e. the integrand of the dy integral if that is the outermost integral */
    FirstOrderArguments args = {{rtols[0], rtols[1]},
                                {atols[0], atols[1]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun};
    args.p = args.Tnu*y;
    args.weight = pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_z_complete, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
};

double integrand_z_complete(double z, FirstOrderArguments args) {
    double Hz = H(z);
    double con = conc(z, args.Mvir_over_Msun);
    args.conc = con;
    args.G0_Gz = Green(0, args.mnu) - Green(z, args.mnu);
    args.Rs = Rs(z, con, Hz, args.Mvir_over_Msun);
    args.Rvir = args.Rs*args.conc;
    args.front_factor = 4.0*std::numbers::pi*_G_*rho0(z, args.Rs, con, args.Mvir_over_Msun)*pow(args.Rs, 3.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight *= args.mnu*args.G0_Gz/(1 + z)/Hz*args.front_factor;
    args.gp = args.p*args.G0_Gz;
    args.r_a = args.rr_here + pow(args.gp, 2.);
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_theta, 0.0, std::numbers::pi, args.rtols[1], args.atols[1], args);
    return I;
}

double zeroth_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu) {
    FirstOrderArguments args = {{rtols[0], rtols[1], rtols[2]},
                                {atols[0], atols[1], atols[2]},
                                N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun};
    
    double I = GaussLaguerre<FirstOrderArguments>(integrand_y_zeroth, args.GaussLaguerreNodes, args);
    return I;
};

double integrand_y_zeroth(double y, FirstOrderArguments args) {
    args.p = args.Tnu*y;
    args.weight = pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_theta_zeroth, 0.0, std::numbers::pi, args.rtols[1], args.atols[1], args);
    return I;
};

double integrand_theta_zeroth(double theta, FirstOrderArguments args) {
    args.weight *= sin(theta);
    return 2.*std::numbers::pi*args.weight;
};

double epsilon_first(double y, double theta, double mass, double z_ini, double rtols[3], double atols[3], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu) {
    FirstOrderArguments args = {{rtols[0], rtols[1], rtols[2]},
                                {atols[0], atols[1], atols[2]},
                                N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun};
    args.p = args.Tnu*y;
    args.theta = theta;
    
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_z_epsilon, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
};

double epsilon_first(double y, double theta, FirstOrderArguments args) {
    args.p = args.Tnu*y;
    args.theta = theta;
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_z_epsilon, 0.0, args.z_ini, args.rtols[0], args.atols[0], args);
    return I;
};

double integrand_z_epsilon(double z, FirstOrderArguments args) {
    args.G0_Gz = Green(0, args.mnu) - Green(z, args.mnu);
    args.gp = args.p*args.G0_Gz;
    args.r_a = args.rr_here + pow(args.gp, 2.);
    args.r = sqrt(args.r_a - 2.*args.gp*args.r_here*cos(args.theta));
    
    double Hz = H(z);
    double con = conc(z, args.Mvir_over_Msun);
    args.Rs = Rs(z, con, Hz, args.Mvir_over_Msun);
    args.Rvir = args.Rs*con;
    args.front_factor = 4.0*std::numbers::pi*_G_*rho0(z, args.Rs, con, args.Mvir_over_Msun)*pow(args.Rs, 3.)*pow(1 + z, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight = args.mnu*args.G0_Gz/Hz/(1 + z);
    
    return args.weight*args.front_factor*LapNFW(args.r, args.Rs);
};
