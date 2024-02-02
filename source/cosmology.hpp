//
//  cosmology.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#ifndef cosmology_hpp
#define cosmology_hpp

#include <stdio.h>
#include <vector>
#include <cmath>
#include <numbers>

#include <iostream>

#define _H0_ 67.66
#define _h_ 0.6766
#define _G_ 6.674e-11
#define _speedoflight_ 2.99792458e+8 // speed of light in m/s
#define _m_to_kpc_ 3.24e-20
#define _H0_kpc_ 2.24688775e-7
#define _OmegaL_ 0.6889
#define _OmegaM_ 0.3111
#define _Msun_ 1.98847e+30 // Solar mass in kg
// Mvir = 3.76e12*1.988e30 // in units of kg
// #define _Mvir_ 4.0365941e+42
// #define _Mvir_over_Msun_ 2.03e+12
#define _sqrt2_ 1.41421356237

double H(double z);
double conc(double z, double Mvir_over_Msun);
double Omega_m(double z);
double rho_crit(double z);
double Rs(double z, double conc, double H, double Mvir_over_Msun);
double rho0(double z, double Rs, double conc, double Mvir_over_Msun);


double Green(double z, double mass);

double LapNFW(double r, double Rs);
double LapNFWKepler(double r, double Rs, double Rvir);

double R1(double r, double Rs);
double R2(double r, double Rs);
double R3(double r, double Rs);



#endif /* cosmology_hpp */
