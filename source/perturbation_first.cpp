//
//  perturbation_first.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include "perturbation_first.hpp"

double first_order(double mass, double z_ini, double rtols[3], double atols[3], double r_here, int N_GaussLaguerre, double Tnu) {
    FirstOrderArguments args = {{rtols[0], rtols[1], rtols[2]},
                                {atols[0], atols[1], atols[2]},
                                N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here};
    
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_z, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
};

double integrand_z(double z, FirstOrderArguments args) {
    double Hz = H(z);
    double con = conc(z);
    args.G0_Gz = Green(0, args.mnu) - Green(z, args.mnu);
    args.Rs = Rs(z, con, Hz);
    args.weight = args.mnu*args.G0_Gz/(1 + z)/Hz*4.0*std::numbers::pi*_G_*rho0(z, args.Rs, con)*pow(args.Rs, 3.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    double I = GaussLaguerre<FirstOrderArguments>(integrand_y, args.GaussLaguerreNodes, args);
    return I;
};

double integrand_y(double y, FirstOrderArguments args) {
    args.p = args.Tnu*y;
    args.gp = args.p*args.G0_Gz;
    args.r_a = args.rr_here + pow(args.gp, 2.);
    args.weight = args.weight*pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrand_theta, 0.0, std::numbers::pi, args.rtols[1], args.atols[1], args);
    return I;
};

double integrand_theta(double theta, FirstOrderArguments args) {
    args.r = sqrt(args.r_a - 2.*args.gp*args.r_here*cos(theta));
    args.weight = args.weight*sin(theta);
    return 2.*std::numbers::pi*args.weight*LapNFW(args.r, args.Rs);
//    auto [I, err] = GaussKronrod(integrand_phi, 0.0, 2*std::numbers::pi, args.rtols[2], args.atols[2], args);
//    return I;
};

/*
double integrand_phi(double phi, FirstOrderArguments args) {
    double G0_Gz_p = args.G0_Gz*args.p;
    double pos_back[3] = {args.x_here[0] - G0_Gz_p*args.sintheta*cos(phi), // Faster to write it like this! Factor ~10 compared to the old
                          args.x_here[1] - G0_Gz_p*args.sintheta*sin(phi),
                          args.x_here[2] - G0_Gz_p*args.costheta};
    return args.phi_integrand;
}
*/
