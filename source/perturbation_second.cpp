//
//  perturbation_second.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 03/03/2023.
//

#include "perturbation_second.hpp"
#include "cosmology.hpp"

double second_order(double mass, double z_ini, double rtols[4], double atols[4], double r_here, int N_GaussLaguerre, double Tnu) {
    SecondOrderArguments args = {{rtols[0], rtols[1], rtols[2], rtols[3]},
                                 {atols[0], atols[1], atols[2], atols[3]},
                                 N_GaussLaguerre, mass, Tnu, r_here, r_here*r_here};
    args.G0 = Green(0, mass);
    args.z_ini = z_ini;
    auto [I, err] = GaussKronrod<SecondOrderArguments>(integrand_z2, 0.0, z_ini, rtols[0], atols[0], args);
    return I;
}

double integrand_z2(double z2, SecondOrderArguments args) {
    double Hz2 = H(z2);
    double conc2 = conc(z2);
    args.Gz2 = Green(z2, args.mnu);
    args.G0_Gz2 = args.G0 - args.Gz2;
    args.Rs2 = Rs(z2, conc2, Hz2);
    args.Rvir2 = args.Rs2*conc2;
    args.front_factor2 = 4.0*std::numbers::pi*_G_*rho0(z2, args.Rs2, conc2)*pow(args.Rs2, 3.)*pow(1 + z2, -2.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
    
    args.weight = args.mnu*args.G0_Gz2/Hz2/(1 + z2);
    auto [I, err] = GaussKronrod<SecondOrderArguments>(integrand_z1, z2, args.z_ini, args.rtols[1], args.atols[1], args);
    return I;
}

double integrand_z1(double z1, SecondOrderArguments args) {
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
    // Future: Use pre-computed GL nodes+weights instead
    double I = GaussLaguerre<SecondOrderArguments>(integrand_y, args.GaussLaguerreNodes, args);
    return I;
}

double integrand_y(double y, SecondOrderArguments args) {
    args.p = args.Tnu*y;
    args.gp1 = args.p*args.G0_Gz1;
    args.gp2 = args.p*args.G0_Gz2;
    args.y_a1 = args.rr_here + pow(args.gp1, 2.);
    args.y_a2 = args.rr_here + pow(args.gp2, 2.);
    args.weight *= pow(args.Tnu, 3.)*pow(y, 2.)/(1.0 + exp(y));
    
    auto [I, err] = GaussKronrod<SecondOrderArguments>(integrand_theta, 0.0, std::numbers::pi, args.rtols[2], args.atols[2], args);
    return I;
}

double integrand_theta(double theta, SecondOrderArguments args) {
    double costheta = cos(theta);
    double y1 = sqrt(args.y_a1 - 2.*args.gp1*args.r_here*costheta);
    double y2 = sqrt(args.y_a2 - 2.*args.gp2*args.r_here*costheta);
    double y1_dot_y2 = args.rr_here + args.gp1*args.gp2 - costheta*args.r_here*(args.gp1 + args.gp2);
    args.weight *= sin(theta);
    double term1 = 0.;
    double term2 = 0.;
    
    // computation of term 1, the (2,2) Laplacian term
    if ((y1 <= args.Rvir1) && (y2 <= args.Rvir2)) {
        // both within Rvir
        // term1 vanishes except when both are within
        term1 = args.G0_Gz1*args.front_factor1*args.front_factor2/(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
    }
    
    // computation of term 2, the (3,1) term
    if (y2 <= args.Rvir2) {
        // term 2 always vanishes if y2 is outside
        if (y1 <= args.Rvir1) {
            // both within Rvir
            term2 = -args.Gz2_Gz1*args.front_factor1*args.front_factor2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*(log((y1 + args.Rs1)/args.Rs1) - y1/(y1 + args.Rs1));
        }
        else {
            // y2 inside, y1 outside; use Kepler for y1
            term2 = -args.Gz2_Gz1*args.front_factor2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*_G_*_Mvir_*pow(1 + args.z1, 3.)/pow(_speedoflight_, 2.)*_m_to_kpc_;
        }
    }
    return 2.*std::numbers::pi*args.weight*(term1 - term2);
}



/*
    CHECKS TO BE MADE
 
    -> Can we reproduce the Kepler potential with the front_factor "formalism"?
    -> -||- for NFW
 
 
 */













/*
 
 double integrand_theta(double theta, SecondOrderArguments args) {
     double costheta = cos(theta);
     double y1 = sqrt(args.y_a1 - 2.*args.gp1*args.r_here*costheta);
     double y2 = sqrt(args.y_a2 - 2.*args.gp2*args.r_here*costheta);
     double y1_dot_y2 = args.rr_here + args.p*(args.p*args.G0_Gz2*args.G0_Gz1 - costheta*args.xGz1pGz2);
     args.weight *= sin(theta);
     
     double R1y1 = args.front_factor1*R1(y1, args.Rs1);
     double R1y2 = args.front_factor2*R1(y2, args.Rs2);
     
     double R2y1 = args.front_factor1*R2(y1, args.Rs1);
     double R2y2 = args.front_factor2*R2(y2, args.Rs2);
     
     double R3y2 = args.front_factor2*R3(y2, args.Rs2);
     
     double term1 = (3.*R1y1 + y1*y1*R2y1)*(3.*R1y2 + y2*y2*R2y2);
     double term2 = -R1y1*(5.*R2y2 + y2*y2*R3y2)*y1_dot_y2;
     
     
     //double front_var = 1./(pow(y1*y2*(y2 + args.Rs2), 3.)*pow(y1 + args.Rs1, 2.));
     double term3 = 0.;
     double term4 = 0.;
     if ((y1 <= args.Rvir1) && (y2 <= args.Rvir2)) {
 //        term3 = pow(y1*y2, 2.)*(y2 + args.Rs2);
         term3 = args.front_factor1*args.front_factor2/(y1*y2*pow((y1 + args.Rs1)*(y2 + args.Rs2), 2.));
     }
     if (y2 <= args.Rvir2) {
         if (y1 <= args.Rvir1) {
             // Both within Rvir
             // term4 = y1_dot_y2*(3*y2 + args.Rs2)*(y1 + args.Rs1)*((y1 + args.Rs1)*log((y1 + args.Rs1)/args.Rs1) - y1);
             term4 = -args.front_factor1*args.front_factor2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*((y1 + args.Rs1)*log((y1 + args.Rs1)/args.Rs1) - y1/(y1 + args.Rs1));
         }
         else {
             // y2 inside, y1 outside; use Kepler for y1
             // term4 = -y1_dot_y2*(3*y2 + args.Rs2)*((y1 + args.Rs1)*log((y1 + args.Rs1)/args.Rs1) - y1);
             
             // Units OK here?
             term4 = -args.front_factor2*y1_dot_y2*(3*y2 + args.Rs2)/pow(y1*y2*(y2 + args.Rs2), 3.)*_G_*_Mvir_*pow(1 + args.z1, 3.)/pow(_speedoflight_, 2.)*_m_to_kpc_;;
         }
     }
     // double all_terms = args.front_factor1*args.front_factor2*front_var*(term3 + term4);
     return std::numbers::pi*args.weight*(term3 - term4);
     //return 2.*std::numbers::pi*args.weight*all_terms;
     // double front_var_reverse = 1./(pow(y2*y1*(y1 + args.Rs1), 3.)*pow(y2 + args.Rs2, 2.));
     // double term3_reverse = pow(y2*y1, 2.)*(y1 + args.Rs1);
     // double term4_reverse = y1_dot_y2*(y2 + args.Rs2)*(3*y1 + args.Rs1)*((y2 + args.Rs2)*log((y2 + args.Rs2)/args.Rs2) - y2);
     // term3 = 0.;
     // all_terms = args.front_factor1*args.front_factor2*front_var_reverse*(term3_reverse + term4_reverse);
     // double Lap1 = pow(y1, 2.)*pow(y1, -3.)*pow(y1 + args.Rs1, 2.);
     // double Lap1_test = LapNFW(y1, args.Rs1);
     // double term3_new = front_var*term3;
     // double term3_test = LapNFW(y1, args.Rs1)*LapNFW(y2, args.Rs2);
     // return std::numbers::pi*args.weight*(term1 + term2);
     //
 }

 */
