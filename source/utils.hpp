//
//  utils.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <vector>

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
#define _Mvir_ 4.0365941e+42
#define _Mvir_over_Msun_ 2.03e+12
#define _sqrt2_ 1.41421356237

std::vector<double> sph2cart(double r, double phi, double sintheta, double costheta);

#endif /* utils_hpp */
