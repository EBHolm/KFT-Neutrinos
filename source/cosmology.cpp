//
//  cosmology.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include "cosmology.hpp"


double H(double z) {
    return _H0_*sqrt(_OmegaM_*pow(1. + z, 3.) + _OmegaL_)/_speedoflight_;
}

double conc(double z, double Mvir_over_Msun) {
    double a = 0.537 + (1.025 - 0.537)*exp(-0.718*pow(z, 1.08));
    double b = -0.097 + 0.024*z;
    return pow(10., a + b*log10(Mvir_over_Msun/(1e+12/_h_)))/0.613;
}

double Omega_m(double z) {
    return _OmegaM_*pow(1. + z, 3.)/(_OmegaM_*pow(1. + z, 3.) + _OmegaL_);
}

double rho_crit(double z) {
    return 3.0*pow(H(z), 2.)/(8.0*std::numbers::pi*_G_*pow(_speedoflight_, -2.)*_m_to_kpc_);
}

double Rs(double z, double conc, double H, double Mvir_over_Msun) {
    double Delta_vir = 18.*pow(std::numbers::pi, 2.) + 82.*(Omega_m(z) - 1.) - 39.*pow(Omega_m(z) - 1., 2.);
    double Rvir = pow(3.0*Mvir_over_Msun*_Msun_/(4.*std::numbers::pi*Delta_vir*rho_crit(z)), 1.0/3.0);
    return Rvir/conc;
}

double rho0(double z, double Rs, double conc, double Mvir_over_Msun) {
    return Mvir_over_Msun*_Msun_/(4.0*std::numbers::pi*pow(Rs, 3.)*pow(1. + z, -3.)*(log(1. + conc) - conc/(1 + conc)));
}

double Green(double z, double mass) {
    /*  Involves evaluating the hypergeometric a=1/2, b=2/3, c=5/3 function approximately on the interval (-40, -5)
        We use a (4,4) MiniMax approximation with (12,12) as machine precision constructed using Horner form */
    double x = -_OmegaM_/_OmegaL_*pow(1. + z, 3.);
    //double num4 = 1. + x*(-0.8484 + x*(0.1531 + (-0.005356 + 0.00001712*x)*x));
    //double denom4 = 1. + x*(-1.048 + x*(0.2695 + (-0.0167 + 0.000165*x)*x));
    double num12 = 0.99999999999999993230472 +
    x*(-3.25057496412430147994808 +x*(4.26290176908055816661688 +
    x*(-2.92054007010084427011892 +x*(1.13616498837807664217507 +
    x*(-0.257142208796449813992234 +x*(0.0336388627043480127418357 +
    x*(-0.00247196836840160076351321 +x*(0.0000969867649764825137762057 +
    x*(-1.87151952467133514384527e-6 +x*(1.54444266660186424724201e-8 +
    (-4.12357970460131970474297e-11 +1.54005371406019611948706e-14*x)*x))))))))));
    double den12 = 1. +x*(-3.4505749641242822609406 +x*(4.85926676190628812550808
    +x*(-3.62572020139819933583841 +x*(1.56274566556175174394472 +
    x*(-0.400038895920748024520823 +x*(0.0606741238537938922181634 +
    x*(-0.00532719506646608896876435 +x*(0.000259422512839402924071888 +
    x*(-6.54781413506950541334025e-6 +x*(7.67661946730659120959631e-8 +
    (-3.42852473695444995156931e-10 +3.67773234079512179756218e-13*x)*x))))))))));
    
    return -pow(1. + z, 2.)/(2.*_H0_kpc_*mass*sqrt(_OmegaL_))*num12/den12;
}

/*
    Potential derivatives
 
 */


double LapNFW(double r, double Rs) {
    return 1.0/(r*pow(r + Rs, 2.));
};

double LapNFWKepler(double r, double Rs, double Rvir) {
    if (r <= Rvir) {
        return 1.0/(r*pow(r + Rs, 2.));
    }
    return 0.;
};

double R1(double r, double Rs) {
    return pow(r, -3.)*(log((r + Rs)/Rs) - r/(r + Rs));
};

double R2(double r, double Rs) {
    return -pow(r, -3.)*(3.*pow(r + Rs, 2.)*log((r + Rs)/Rs) - (r*(4.*r + 3.*Rs)))/pow(r*(r + Rs), 2.);
};

double R3(double r, double Rs) {
    return pow(r, -3.)*(15.*pow(r + Rs, 3.)*log((r + Rs)/Rs) - r*(23.*r*r + 36.*r*Rs + 15.*Rs*Rs))/(pow(r, 4.)*pow(r + Rs, 3.));
};
