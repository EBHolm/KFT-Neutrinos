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
#include "utils.hpp"

double H(double z);
double conc(double z);
double Omega_m(double z);
double rho_crit(double z);
double Rs(double z, double conc, double H);
double rho0(double z, double Rs, double conc);


double Green(double z, double mass);

double LapNFW(double r, double Rs);
double LapNFWKepler(double r, double Rs, double Rvir);

double R1(double r, double Rs);
double R2(double r, double Rs);
double R3(double r, double Rs);



#endif /* cosmology_hpp */
