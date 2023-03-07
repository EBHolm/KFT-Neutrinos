//
//  perturbation_second.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 03/03/2023.
//

#include "perturbation_second.hpp"

double second_order(double mass, double z_ini, double rtols[4], double atols[4], double r_here, int N_GaussLaguerre, double Tnu) {
    SecondOrderArguments args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here};
    // TO-DO: Set args
    auto [I, err] = GaussKronrod(integrand_z2, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
}

double integrand_z2(double z2, SecondOrderArguments args) {
    // TO-DO: Set new args
    // Make sure that this integral is actually from 0 to z2, not from 3 to z2!
    auto [I, err] = GaussKronrod(integrand_z1, 0.0, z2, args.rtols[1], args.atols[1], args);
    return I;
}

double integrand_z1(double z1, SecondOrderArguments args) {
    // TO-DO: Set new args
    // Future: Use pre-computed GL nodes+weights instead
    double I = GaussLaguerre(integrand_y, args.GaussLaguerreNodes, args);
    return I;
}

double integrand_y(double y, SecondOrderArguments args) {
    // TO-DO: Set new args
    auto [I, err] = GaussKronrod(integrand_theta, 0.0, std::numbers::pi, args.rtols[2], args.atols[2], args);
    return I;
}

double integrand_theta(double theta, SecondOrderArguments args) {
    // TO-DO: Set new args
    //        Check analytically whether explicit phi integrand is needed
    auto [I, err] = GaussKronrod(integrand_phi, 0.0, 2.*std::numbers::pi, args.rtols[3], args.atols[3], args);
    return I;
}

double integrand_phi(double phi, SecondOrderArguments args) {
    return 1.0;
}
