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
#include "perturbation_second.hpp"

int main(int argc, const char * argv[]) {
    /*
    if (argc < 2) {
        throw std::invalid_argument("You must give the mass as command line argument.");
    }
    const double mass = std::stod(argv[1]);
    */
    
    double mass = 0.1;
    double r_here = 8.2;
    
    double rtols[3] = {1e-6, 1e-4, 1e-4};
    double atols[3] = {1e-35, 1e-35, 1e-35};
    
    double rtols_2[4] = {1e-6, 1e-4, 1e-4, 1e-4};
    double atols_2[4] = {1e-35, 1e-35, 1e-35, 1e-35};
    
    int GaussLaguerreNodes = 50;
    
    double analytical_free = 1.06738178968e-10;
    double first = first_order(mass, 3.0, rtols, atols, r_here, GaussLaguerreNodes);
    
    /* Should give 1.2155 at m=0.1 eV */
    /*         and 1.0018 at m=0.01 eV */
    std::cout << "First order: " << first << ", corresponding to a clustering factor " << 1+first/analytical_free << " with m=" << mass << " eV.\n";
    
    
    double second = second_order(mass, 3.0, rtols_2, atols_2, r_here, GaussLaguerreNodes);
    std::cout << "Second order: " << second << ", corresponding to a clustering factor " << 1+(first + second)/analytical_free << " with m=" << mass << " eV.\n";
    
    return 0;
}
