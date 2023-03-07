//
//  utils.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include "utils.hpp"
#include <vector>
#include <cmath>

std::vector<double> sph2cart(double r, double phi, double sintheta, double costheta) {
    double x1 = r*sintheta*cos(phi);
    double x2 = r*sintheta*sin(phi);
    double x3 = r*costheta;
    return {x1, x2, x3};
};
