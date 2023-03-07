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
    if (argc < 2) {
        throw std::invalid_argument("You must give the mass as command line argument.");
    }
    
    const double mass = std::stod(argv[1]);
    
    // double mass = 0.1;
    
    // std::cout << "Computing first order perturbation theory result with mass = " << mass << " eV.\n";
    
    double rtols[3] = {1e-6, 1e-4, 1e-4};
    double atols[3] = {1e-35, 1e-35, 1e-35};
    double rtols_2[4] = {1e-6, 1e-4, 1e-4, 1e-4};
    double atols_2[4] = {1e-35, 1e-35, 1e-35, 1e-35};
    double x_here[3] = {0.0, 0.0, 8.0}; // For optimization purposes, the code assumes that x_here only has a z-component
    int GaussLaguerreNodes = 50;
    
    double R = Rs(0.0, conc(0.0), H(0.0));
    double asdf = rho0(0.0, 19.9, conc(0.0));
    
    /*
    std::function<double(double, FunctionArguments)> test_phi = integrand_phi;
    FunctionArguments args = {
        .rtols = {1e-2, 1e-2, 1e-2},   // for z, theta and phi integrals
        .atols = {1e-25, 1e-25, 1e-25},
        .GaussLaguerreNodes = 20,      // for p integral
        
        .x_here = {8.0, 0.0, 0.0},
        .mnu = mass,                   // eV
        .Tnu = 0.0001676375864435959,  // eV
        
        .G0_Gz = 1.0,
        .weight = 1.0,
        .Rs = 1.0,
        .p = 1.0,
    };
    // std::cout << "Call to integrand_phi = " << test_phi(0.5, args) << "\n";
    std::function<double(double, FunctionArguments)> FermiDirac = [](double x, FunctionArguments args) {return 1.0/(1.0 + exp(x));};
    double FermiDirac_GL = GaussLaguerre(FermiDirac, GaussLaguerreNodes, args); // 0.69315, correct!
    auto [FermiDirac_GK, err] = GaussKronrod(FermiDirac, 0, 20, 1e-12, 1e-35, args);
    std::cout << "FermiDirac integral: GK=" << FermiDirac_GK << ", GL=" << FermiDirac_GL << "\n";*/
    
    double analytical_free = 1.06738178968e-10;
    double first = first_order(mass, 3.0, rtols, atols, 8.0, GaussLaguerreNodes);
    
    /* Should give 1.2155 at m=0.1 eV */
    /*         and 1.0018 at m=0.01 eV */
    std::cout << "First order: " << first << ", corresponding to a clustering factor " << 1+first/analytical_free << " with m=" << mass << " eV.\n";
    
    
    double second = second_order(mass, 3.0, rtols_2, atols_2, 8.0, GaussLaguerreNodes);
    std::cout << "Second order: " << second << ", corresponding to a clustering factor " << 1+(first + second)/analytical_free << " with m=" << mass << " eV.\n";
    
    return 0;
}
