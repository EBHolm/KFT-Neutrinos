//
//  perturbation_second_vlasov.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 03/03/2023.
//

#include "perturbation_second_vlasov.hpp"
#include "cosmology.hpp"

double second_order_vlasov(double mass, double z_ini, double rtols[4], double atols[4], double r_here, int N_GaussLaguerre, double Tnu) {
    SecondOrderArgumentsVlasov args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here};
    args.G0 = Green(0, mass);
    args.z_ini = z_ini;
    auto [I, err] = GaussKronrod<SecondOrderArgumentsVlasov>(integrand_z2_vlasov, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
}

double integrand_z2_vlasov(double z2, SecondOrderArgumentsVlasov args) {
    double Hz2 = H(z2);
    double conc2 = conc(z2);
    args.Gz2 = Green(z2, args.mnu);
    args.G0_Gz2 = args.G0 - args.Gz2;
    args.Rs2 = Rs(z2, conc2, Hz2);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(z2, args.Rs2, conc2)*pow(args.Rs2, 3.)*pow(1 + z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    
    args.weight = args.mnu*args.G0_Gz2/Hz2/(1 + z2);
    auto [I, err] = GaussKronrod<SecondOrderArgumentsVlasov>(integrand_z1_vlasov, z2, args.z_ini, args.rtols[1], args.atols[1], args);
    return I;
}

double integrand_z1_vlasov(double z1, SecondOrderArgumentsVlasov args) {
    double Hz1 = H(z1);
    double conc1 = conc(z1);
    double Gz1 = Green(z1, args.mnu);
    args.z1 = z1;
    args.G0_Gz1 = args.G0 - Gz1;
    args.Gz2_Gz1 = args.Gz2 - Gz1;
    args.Rs1 = Rs(z1, conc1, Hz1);
    args.Rvir1 = args.Rs1*conc1;
    args.front_factor1 = 4.0*std::numbers::pi*_G_*rho0(z1, args.Rs1, conc1)*pow(args.Rs1, 3.)*pow(1 + z1, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;

    args.weight *= args.mnu*args.G0_Gz1/Hz1/(1 + z1);
    // Future: Use pre-computed GL nodes+weights instead
    double I = GaussLaguerre<SecondOrderArgumentsVlasov>(integrand_y_vlasov, args.GaussLaguerreNodes, args);
    return I;
}

double integrand_y_vlasov(double y, SecondOrderArgumentsVlasov args) {
    args.p = args.Tnu*y;
    args.gp1 = args.p*args.G0_Gz1;
    args.gp2 = args.p*args.G0_Gz2;
    args.y_a1 = args.rr_here + pow(args.gp1, 2.);
    args.y_a2 = args.rr_here + pow(args.gp2, 2.);
    args.weight *= pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    
    auto [I, err] = GaussKronrod<SecondOrderArgumentsVlasov>(integrand_theta_vlasov, 0.0, std::numbers::pi, args.rtols[2], args.atols[2], args);
    return I;
}

double integrand_theta_vlasov(double theta, SecondOrderArgumentsVlasov args) {
    double costheta = cos(theta);
    double y1 = sqrt(args.y_a1 - 2.*args.gp1*args.r_here*costheta);
    double y2 = sqrt(args.y_a2 - 2.*args.gp2*args.r_here*costheta);
    double y1_dot_y2 = args.rr_here + args.gp1*args.gp2 - costheta*args.r_here*(args.gp1 + args.gp2);
    args.weight *= sin(theta);
    
    // computation of term 1, the (2,2) Laplacian term
    double term1 = 1./(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
    
    // computation of term 2, the (3,1) term
    double term2 = -y1_dot_y2*(3*y1 + args.Rs1)/pow(y2*y1*(y1 + args.Rs1), 3.)*(log((y2 + args.Rs2)/args.Rs2) - y2/(y2 + args.Rs2));
    
    return 2.*std::numbers::pi*args.weight*args.front_factor1*args.front_factor2*(term1 - term2);
}

double integrand_z2z1_vlasov(double z2, double z1, double mass, double z_ini, double rtols[4], double atols[4], double r_here, int N_GaussLaguerre, double Tnu) {
    /* Used for plotting the integrand in the (z2, z1)-plane */
    
    // General args
    SecondOrderArgumentsVlasov args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here};
    args.G0 = Green(0, mass);
    args.z_ini = z_ini;
    
    // z2 args
    double Hz2 = H(z2);
    double conc2 = conc(z2);
    args.Gz2 = Green(z2, args.mnu);
    args.G0_Gz2 = args.G0 - args.Gz2;
    args.Rs2 = Rs(z2, conc2, Hz2);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(z2, args.Rs2, conc2)*pow(args.Rs2, 3.)*pow(1 + z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight = args.mnu*args.G0_Gz2/Hz2/(1 + z2);
    
    // z1 args
    double Hz1 = H(z1);
    double conc1 = conc(z1);
    double Gz1 = Green(z1, args.mnu);
    args.z1 = z1;
    args.G0_Gz1 = args.G0 - Gz1;
    args.Gz2_Gz1 = args.Gz2 - Gz1;
    args.Rs1 = Rs(z1, conc1, Hz1);
    args.Rvir1 = args.Rs1*conc1;
    args.front_factor1 = 4.0*std::numbers::pi*_G_*rho0(z1, args.Rs1, conc1)*pow(args.Rs1, 3.)*pow(1 + z1, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight *= args.mnu/Hz1/(1 + z1);
    
    double I = GaussLaguerre<SecondOrderArgumentsVlasov>(integrand_y_vlasov, args.GaussLaguerreNodes, args);
    return I;
}

