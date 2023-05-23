//
//  main.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include <iostream>
#include <stdexcept>
#include <string>
#include "perturbation_first.hpp"

int main(int argc, const char * argv[]) {
    double mass = 0.01;
    
    std::cout << "Computing first order perturbation theory result with mass = " << mass << " eV.\n";
    
    double rtols[3] = {1e-7, 1e-5, 1e-5};
    double atols[3] = {1e-35, 1e-35, 1e-35};
    int GaussLaguerreNodes = 70;
    
    double analytical_free = 1.06738178968e-10;
    double first = first_order(mass, 3.0, rtols, atols, 8.0, GaussLaguerreNodes);
    
    std::cout << "First order: " << first << ", corresponding to a clustering factor " << 1+first/analytical_free << " with m=" << mass << " eV.\n";
    return 0;
}
