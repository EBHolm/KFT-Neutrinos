//
//  perturbation_second_kft.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 03/03/2023.
//

#include "perturbation_second_kft.hpp"
#include "cosmology.hpp"

double second_order_kft(double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag, double Tnu) {
    SecondOrderArgumentsKFT args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun, terms_flag};
    args.G0 = Green(0, mass);
    args.z_ini = z_ini;
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_z2_kft, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
}

double integrand_z2_kft(double z2, SecondOrderArgumentsKFT args) {
    args.z2 = z2;
    double Hz2 = H(z2);
    double conc2 = conc(z2, args.Mvir_over_Msun);
    args.Gz2 = Green(z2, args.mnu);
    args.G0_Gz2 = args.G0 - args.Gz2;
    args.Rs2 = Rs(z2, conc2, Hz2, args.Mvir_over_Msun);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(z2, args.Rs2, conc2, args.Mvir_over_Msun)*pow(args.Rs2, 3.)*pow(1 + z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight = args.mnu*args.G0_Gz2/Hz2/(1 + z2)*args.front_factor2;
    
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_z1_kft, z2, args.z_ini, args.rtols[1], args.atols[1], args);
    // auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_y_kft, 1e-8, 20.0, args.rtols[2], args.atols[2], args);
    return I;
}

double integrand_z1_kft(double z1, SecondOrderArgumentsKFT args) {
    double Hz1 = H(z1);
    double conc1 = conc(z1, args.Mvir_over_Msun);
    double Gz1 = Green(z1, args.mnu);
    args.z1 = z1;
    args.G0_Gz1 = args.G0 - Gz1;
    args.Rs1 = Rs(z1, conc1, Hz1, args.Mvir_over_Msun);
    args.Rvir1 = args.Rs1*conc1;
    args.front_factor1 = 4.0*std::numbers::pi*_G_*rho0(z1, args.Rs1, conc1, args.Mvir_over_Msun)*pow(args.Rs1, 3.)*pow(1 + z1, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;

    args.weight *= args.mnu/Hz1/(1 + z1)*args.front_factor1;
    // Future: Use pre-computed GL nodes+weights instead
    double I = GaussLaguerre<SecondOrderArgumentsKFT>(integrand_y_kft, args.GaussLaguerreNodes, args);
    //auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_y_kft, 1e-8, 20.0, args.rtols[2], args.atols[2], args);
    return I;
}

double integrand_y_kft(double y, SecondOrderArgumentsKFT args) {
    args.p = args.Tnu*y;
    args.gp2 = args.p*args.G0_Gz2;
    args.y_a2 = args.rr_here + pow(args.gp2, 2.);
    args.gp1 = args.p*args.G0_Gz1;
    args.y_a1 = args.rr_here + pow(args.gp1, 2.);
    args.weight *= pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_theta_kft, 0.0, std::numbers::pi, args.rtols[2], args.atols[2], args);
    // auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_theta_kft_term3, 0.0, std::numbers::pi, args.rtols[2], args.atols[2], args);
    return I;
}

double integrand_theta_kft(double theta, SecondOrderArgumentsKFT args) {
    double costheta = cos(theta);
    double y1 = sqrt(args.y_a1 - 2.*args.gp1*args.r_here*costheta);
    double y2 = sqrt(args.y_a2 - 2.*args.gp2*args.r_here*costheta);
    double y1_dot_y2 = args.rr_here + args.gp1*args.gp2 - costheta*args.r_here*(args.gp1 + args.gp2);
    args.weight *= sin(theta);
    
    // Only computing NFW part for now
    // This is (A.19) in the notes (as of dec. 12)
    // But note the replacement z1 <-> z2 between here and (A.19)
    
    double term1 = 0.;
    double term2 = 0.;
    double term3 = 0.;
    double term4 = 0.;
    
    if (args.terms_flag == 0) {
        // Compute all terms
        // computation of term 1, the (3,1) term
        term1 = -args.G0_Gz1*y1_dot_y2*(3*y1 + args.Rs1)/pow(y2*y1*(y1 + args.Rs1), 3.)*(log((y2 + args.Rs2)/args.Rs2) - y2/(y2 + args.Rs2));

        // computation of term 2, the (1,3) term
        term2 = -args.G0_Gz2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*(log((y1 + args.Rs1)/args.Rs1) - y1/(y1 + args.Rs1));
        
        // computation of term 3, the (2,2) Laplacian term
        term3 = args.G0_Gz1/(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
        
        // computation of term 4, the (2,2) cross-contracted term
        double Log1 = log(1 + y1/args.Rs1);
        double Log2 = log(1 + y2/args.Rs2);
        double cross = y1*args.Rs1*(2*Log1 - 1)*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*y1_dot_y2*(3*Log2 - 4) +                                               // first line done
                                    3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2)) +              // second line done
                       pow(args.Rs1, 2.)*Log1*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*y1_dot_y2*(3*Log2 - 4) +                                                 // third line
                                                3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2)) +  // fourth line, exactly same as second line
                       y1*y1*( pow(y1, 2.)*pow(y2, 4.)*((1 - 2*Log2) + Log1*(3*Log2 - 2)) +                                                                             // fifth line
                               pow(y2, 2.)*y1_dot_y2*y1_dot_y2*(3*Log1 - 4)*(3*Log2 - 4) +                                                                                        // sixth line
                               y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.)*(3*Log1 - 2) + 3*y1_dot_y2*y1_dot_y2*(3*Log1 - 4)) +                                                      // seventh line
                               pow(args.Rs2, 2.)*Log2*(pow(y1*y2, 2.)*(3*Log1 - 2)  + 3*y1_dot_y2*y1_dot_y2*(3*Log1 - 4)));                                                       // eigth line
        double denom = pow(y1*y2, 5.)*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.);

        term4 = args.G0_Gz2*cross/denom;
    }
    else if (args.terms_flag == 1) {
        term1 = -args.G0_Gz1*y1_dot_y2*(3*y1 + args.Rs1)/pow(y2*y1*(y1 + args.Rs1), 3.)*(log((y2 + args.Rs2)/args.Rs2) - y2/(y2 + args.Rs2));
    }
    else if (args.terms_flag == 2) {
        term2 = -args.G0_Gz2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*(log((y1 + args.Rs1)/args.Rs1) - y1/(y1 + args.Rs1));
    }
    else if (args.terms_flag == 3) {
        term3 = args.G0_Gz1/(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
    }
    else if (args.terms_flag == 4) {
        double Log1 = log(1 + y1/args.Rs1);
        double Log2 = log(1 + y2/args.Rs2);
        double cross = y1*args.Rs1*(2*Log1 - 1)*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*(3*Log2 - 4) +                                               // first line done
                                    3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2)) +              // second line done
                       pow(args.Rs1, 2.)*Log1*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*(3*Log2 - 4) +                                                 // third line
                                                3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2)) +  // fourth line, exactly same as second line
                       y1*y1*( pow(y1, 2.)*pow(y2, 4.)*((1 - 2*Log2) + Log1*(3*Log2 - 2)) +                                                                             // fifth line
                               pow(y2, 2.)*y1_dot_y2*(3*Log1 - 4)*(4*Log2 - 4) +                                                                                        // sixth line
                               y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.)*(3*Log1 - 2) + 3*y1_dot_y2*(3*Log1 - 4)) +                                                      // seventh line
                               pow(args.Rs2, 2.)*Log2*(pow(y1*y2, 2.)*(3*Log1 - 2)  + 3*y1_dot_y2*(3*Log1 - 4)));                                                       // eigth line
        double denom = pow(y1*y2, 5.)*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.);

        term4 = args.G0_Gz2*cross/denom;
    }
    
    return 2.*std::numbers::pi*args.weight*(term1 + term2 + term3 + term4);
}

double integrand_theta_kft_term3(double theta, SecondOrderArgumentsKFT args) {
    //args.weight *= sin(theta);
    args.theta = theta;
    double Hz2 = H(args.z2);
    double conc2 = conc(args.z2, args.Mvir_over_Msun);
    double Gz2 = Green(args.z2, args.mnu);
    args.G0_Gz2 = args.G0 - Gz2;
    args.Rs2 = Rs(args.z2, conc2, Hz2, args.Mvir_over_Msun);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(args.z2, args.Rs2, conc2, args.Mvir_over_Msun)*pow(args.Rs2, 3.)*pow(1 + args.z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.gp2 = args.p*args.G0_Gz2;
    args.y_a2 = args.rr_here + pow(args.gp2, 2.);
    double y2 = sqrt(args.y_a2 - 2.*args.gp2*args.r_here*cos(args.theta));
    args.integrand_2 = 0.5*args.mnu*args.G0_Gz2/Hz2/(1 + args.z2)*args.front_factor2*LapNFW(y2, args.Rs2);
    
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_z1_kft_term3, 0.0, args.z_ini, args.rtols[1], args.atols[1], args);
    return I;
};

double integrand_z1_kft_term3(double z1, SecondOrderArgumentsKFT args) {
    double Hz1 = H(z1);
    double conc1 = conc(z1, args.Mvir_over_Msun);
    double Gz1 = Green(z1, args.mnu);
    args.z1 = z1;
    args.G0_Gz1 = args.G0 - Gz1;
    args.Rs1 = Rs(z1, conc1, Hz1, args.Mvir_over_Msun);
    args.Rvir1 = args.Rs1*conc1;
    args.front_factor1 = 4.0*std::numbers::pi*_G_*rho0(z1, args.Rs1, conc1, args.Mvir_over_Msun)*pow(args.Rs1, 3.)*pow(1 + z1, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.gp1 = args.p*args.G0_Gz1;
    args.y_a1 = args.rr_here + pow(args.gp1, 2.);
    double y1 = sqrt(args.y_a1 - 2.*args.gp1*args.r_here*cos(args.theta));
    args.integrand_1 = args.mnu/Hz1/(1 + z1)*args.G0_Gz1*args.front_factor1*LapNFW(y1, args.Rs1);
    // printf("integrand1 = %g, integrand2 = %g\n", args.integrand_1, args.integrand_2);
    
    double y = args.p/args.Tnu;
    double y_weight = pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    double angular_weight = 2.*std::numbers::pi*sin(args.theta);
    
    return angular_weight*y_weight*args.integrand_2*args.integrand_1;
}

double integrand_z2z1_kft(double z2, double z1, double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, double Tnu) {
    /* Used for plotting the integrand in the (z2, z1)-plane */
    
    // General args
    SecondOrderArgumentsKFT args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun};
    args.G0 = Green(0, mass);
    args.z_ini = z_ini;
    
    // z2 args
    double Hz2 = H(z2);
    double conc2 = conc(z2, args.Mvir_over_Msun);
    args.Gz2 = Green(z2, args.mnu);
    args.G0_Gz2 = args.G0 - args.Gz2;
    args.Rs2 = Rs(z2, conc2, Hz2, args.Mvir_over_Msun);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(z2, args.Rs2, conc2, args.Mvir_over_Msun)*pow(args.Rs2, 3.)*pow(1 + z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight = args.mnu*args.G0_Gz2/Hz2/(1 + z2);
    
    // z1 args
    double Hz1 = H(z1);
    double conc1 = conc(z1, args.Mvir_over_Msun);
    double Gz1 = Green(z1, args.mnu);
    args.z1 = z1;
    args.G0_Gz1 = args.G0 - Gz1;
    args.Rs1 = Rs(z1, conc1, Hz1, args.Mvir_over_Msun);
    args.Rvir1 = args.Rs1*conc1;
    args.front_factor1 = 4.0*std::numbers::pi*_G_*rho0(z1, args.Rs1, conc1, args.Mvir_over_Msun)*pow(args.Rs1, 3.)*pow(1 + z1, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.weight *= args.mnu/Hz1/(1 + z1);
    
    double I = GaussLaguerre<SecondOrderArgumentsKFT>(integrand_y_kft, args.GaussLaguerreNodes, args);
    return I;
}


double integrate_epsilon_first(double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag, double Tnu) {
    FirstOrderArguments args = {{rtols[0], rtols[1], rtols[2]},
                                {atols[0], atols[1], atols[2]},
                                N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun};
    args.z_ini = z_ini;
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrate_epsilon_y_first, 1e-8, 15.0, rtols[0], atols[0], args);
    return I;
};

double integrate_epsilon_y_first(double y, FirstOrderArguments args) {
    args.p = y*args.Tnu;
    auto [I, err] = GaussKronrod<FirstOrderArguments>(integrate_epsilon_theta_first, 0.0, std::numbers::pi, args.rtols[0], args.atols[0], args);
    return I;
};

double integrate_epsilon_theta_first(double theta, FirstOrderArguments args) {
    args.theta = theta;
    double y = args.p/args.Tnu;
    double y_weight = pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    double angular_weight = 2.*std::numbers::pi*sin(args.theta);
    
    double eps_first = epsilon_first(args.p/args.Tnu, args.theta, args);
    return angular_weight*y_weight*0.5*eps_first*eps_first;
};


double integrate_epsilon_second(double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag, double Tnu) {
    SecondOrderArgumentsKFT args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun, terms_flag};
    args.z_ini = z_ini;
    args.G0 = Green(0., args.mnu);
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrate_epsilon_y_second, 1e-8, 25.0, rtols[0], atols[0], args);
    return I;
};

double integrate_epsilon_y_second(double y, SecondOrderArgumentsKFT args) {
    args.p = y*args.Tnu;
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrate_epsilon_theta_second, 0.0, std::numbers::pi, args.rtols[0], args.atols[0], args);
    return I;
};

double integrate_epsilon_theta_second(double theta, SecondOrderArgumentsKFT args) {
    args.theta = theta;
    double y = args.p/args.Tnu;
    double y_weight = pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    double angular_weight = 2.*std::numbers::pi*sin(args.theta);
    
    double eps_second = epsilon_second(args.p/args.Tnu, args.theta, args);
    return angular_weight*y_weight*eps_second;
};

double epsilon_second(double y, double theta, double mass, double z_ini, double rtols[4], double atols[4], double r_here, double Mvir_over_Msun, int N_GaussLaguerre, int terms_flag, double Tnu) {
    SecondOrderArgumentsKFT args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here, Mvir_over_Msun, terms_flag};
    args.z_ini = z_ini;
    args.p = args.Tnu*y;
    args.theta = theta;
    args.G0 = Green(0., args.mnu);
    
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_z2_epsilon, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
}

double epsilon_second(double y, double theta,  SecondOrderArgumentsKFT args) {
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_z2_epsilon, 0.0, args.z_ini, args.rtols[0], args.atols[0], args);
    return I;
}

double integrand_z2_epsilon(double z2, SecondOrderArgumentsKFT args) {
    double Hz2 = H(z2);
    double conc2 = conc(z2, args.Mvir_over_Msun);
    double Gz2 = Green(z2, args.mnu);
    args.G0_Gz2 = args.G0 - Gz2;
    args.Rs2 = Rs(z2, conc2, Hz2, args.Mvir_over_Msun);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(z2, args.Rs2, conc2, args.Mvir_over_Msun)*pow(args.Rs2, 3.)*pow(1 + z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.gp2 = args.p*args.G0_Gz2;
    args.y_a2 = args.rr_here + pow(args.gp2, 2.);
    args.weight = args.mnu*args.G0_Gz2/Hz2/(1 + z2)*args.front_factor2;
    
    auto [I, err] = GaussKronrod<SecondOrderArgumentsKFT>(integrand_z1_epsilon, z2, args.z_ini, args.rtols[0], args.atols[0], args);
    return I;
};

double integrand_z1_epsilon(double z1, SecondOrderArgumentsKFT args) {
    double Hz1 = H(z1);
    double conc1 = conc(z1, args.Mvir_over_Msun);
    double Gz1 = Green(z1, args.mnu);
    args.z1 = z1;
    args.G0_Gz1 = args.G0 - Gz1;
    args.Rs1 = Rs(z1, conc1, Hz1, args.Mvir_over_Msun);
    args.Rvir1 = args.Rs1*conc1;
    args.front_factor1 = 4.0*std::numbers::pi*_G_*rho0(z1, args.Rs1, conc1, args.Mvir_over_Msun)*pow(args.Rs1, 3.)*pow(1 + z1, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    args.gp1 = args.p*args.G0_Gz1;
    args.y_a1 = args.rr_here + pow(args.gp1, 2.);
    args.weight *= args.mnu/Hz1/(1 + z1)*args.front_factor1;
    
    double costheta = cos(args.theta);
    double y1 = sqrt(args.y_a1 - 2.*args.gp1*args.r_here*costheta);
    double y2 = sqrt(args.y_a2 - 2.*args.gp2*args.r_here*costheta);
    double y1_dot_y2 = args.rr_here + args.gp1*args.gp2 - costheta*args.r_here*(args.gp1 + args.gp2);
    
    double term1 = 0.;
    double term2 = 0.;
    double term3 = 0.;
    double term4 = 0.;
    
    if (args.terms_flag == 0) {
        // Compute all terms
        // computation of term 1, the (3,1) term
        term1 = -args.G0_Gz1*y1_dot_y2*(3*y1 + args.Rs1)/pow(y2*y1*(y1 + args.Rs1), 3.)*(log((y2 + args.Rs2)/args.Rs2) - y2/(y2 + args.Rs2));

        // computation of term 2, the (1,3) term
        term2 = -args.G0_Gz2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*(log((y1 + args.Rs1)/args.Rs1) - y1/(y1 + args.Rs1));
        
        // computation of term 3, the (2,2) Laplacian term
        term3 = args.G0_Gz1/(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
        
        // computation of term 4, the (2,2) cross-contracted term
        double Log1 = log(1 + y1/args.Rs1);
        double Log2 = log(1 + y2/args.Rs2);
        double cross = y1*args.Rs1*(2*Log1 - 1)*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*y1_dot_y2*(3*Log2 - 4) +                                               // first line done
                                    3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2)) +              // second line done
                       pow(args.Rs1, 2.)*Log1*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*y1_dot_y2*(3*Log2 - 4) +                                                 // third line
                                                3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2*y1_dot_y2)) +  // fourth line, exactly same as second line
                       y1*y1*( pow(y1, 2.)*pow(y2, 4.)*((1 - 2*Log2) + Log1*(3*Log2 - 2)) +                                                                             // fifth line
                               pow(y2, 2.)*y1_dot_y2*y1_dot_y2*(3*Log1 - 4)*(3*Log2 - 4) +                                                                                        // sixth line
                               y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.)*(3*Log1 - 2) + 3*y1_dot_y2*y1_dot_y2*(3*Log1 - 4)) +                                                      // seventh line
                               pow(args.Rs2, 2.)*Log2*(pow(y1*y2, 2.)*(3*Log1 - 2)  + 3*y1_dot_y2*y1_dot_y2*(3*Log1 - 4)));                                                       // eigth line
        double denom = pow(y1*y2, 5.)*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.);

        term4 = args.G0_Gz2*cross/denom;
    }
    else if (args.terms_flag == 1) {
        term1 = -args.G0_Gz1*y1_dot_y2*(3*y1 + args.Rs1)/pow(y2*y1*(y1 + args.Rs1), 3.)*(log((y2 + args.Rs2)/args.Rs2) - y2/(y2 + args.Rs2));
    }
    else if (args.terms_flag == 2) {
        term2 = -args.G0_Gz2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*(log((y1 + args.Rs1)/args.Rs1) - y1/(y1 + args.Rs1));
    }
    else if (args.terms_flag == 3) {
        term3 = args.G0_Gz1/(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
    }
    else if (args.terms_flag == 4) {
        double Log1 = log(1 + y1/args.Rs1);
        double Log2 = log(1 + y2/args.Rs2);
        double cross = y1*args.Rs1*(2*Log1 - 1)*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*(3*Log2 - 4) +                                               // first line done
                                    3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2)) +              // second line done
                       pow(args.Rs1, 2.)*Log1*( pow(y1, 2.)*pow(y2, 4.)*(3*Log2 - 2) + 3*y2*y2*y1_dot_y2*(3*Log2 - 4) +                                                 // third line
                                                3*y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.) + 3*y1_dot_y2) + 3*Log2*pow(args.Rs2, 2.)*(pow(y1*y2, 2.) + 3*y1_dot_y2)) +  // fourth line, exactly same as second line
                       y1*y1*( pow(y1, 2.)*pow(y2, 4.)*((1 - 2*Log2) + Log1*(3*Log2 - 2)) +                                                                             // fifth line
                               pow(y2, 2.)*y1_dot_y2*(3*Log1 - 4)*(4*Log2 - 4) +                                                                                        // sixth line
                               y2*args.Rs2*(2*Log2 - 1)*(pow(y1*y2, 2.)*(3*Log1 - 2) + 3*y1_dot_y2*(3*Log1 - 4)) +                                                      // seventh line
                               pow(args.Rs2, 2.)*Log2*(pow(y1*y2, 2.)*(3*Log1 - 2)  + 3*y1_dot_y2*(3*Log1 - 4)));                                                       // eigth line
        double denom = pow(y1*y2, 5.)*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.);

        term4 = args.G0_Gz2*cross/denom;
    }
    
    return args.weight*(term1 + term2 + term3 + term4);
};
