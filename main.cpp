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
#include "perturbation_second_kft.hpp"

int main(int argc, const char * argv[]) {
    /*
    if (argc < 2) {
        throw std::invalid_argument("You must give the mass as command line argument.");
    }
    const double mass = std::stod(argv[1]);
    */
    
    
    double mass = 0.1;
    double r_here = 8.2;
    
    double rtols[3] = {1e-3, 1e-3, 1e-3};
    double rtols_2[4] = {1e-3, 1e-3, 1e-3, 1e-3};
    
    double atols[3] = {1e-35, 1e-35, 1e-35};
    double atols_2[4] = {1e-35, 1e-35, 1e-35, 1e-35};
    
    double Mvir_over_Msun = 2.03e+12;
    
    int GaussLaguerreNodes = 50;
    
    double analytical_free = 1.06738178968e-10;
    double zeroth = zeroth_order(mass, 3.0, rtols, atols, r_here, Mvir_over_Msun, GaussLaguerreNodes);
    double first = first_order(mass, 3.0, rtols, atols, r_here, Mvir_over_Msun, GaussLaguerreNodes);
    std::cout << "At Mvir = " << Mvir_over_Msun << " Msun:\n";
    std::cout << "Zeroth order: " << zeroth << ", corresponding to a clustering factor " << zeroth/analytical_free << " with m=" << mass << " eV.\n";
    std::cout << "First order: " << first << ", corresponding to a clustering factor " << 1+first/analytical_free << " with m=" << mass << " eV.\n";
    
    double second_kft = second_order_kft(mass, 3.0, rtols_2, atols_2, r_here, Mvir_over_Msun, GaussLaguerreNodes);
    std::cout << "Second order KFT: " << second_kft << ", corresponding to a clustering factor " << 1+(first + second_kft)/analytical_free << " with m=" << mass << " eV.\n\n";
    
    double integrated_first = integrate_epsilon_first(mass, 3.0, rtols, atols, r_here, Mvir_over_Msun, GaussLaguerreNodes);
    double integrated_second = integrate_epsilon_second(mass, 3.0, rtols_2, atols_2, r_here, Mvir_over_Msun, GaussLaguerreNodes, 3);
    std::cout << "\nI from first:  " << integrated_first/analytical_free << "\n";
    std::cout << "0.5*I^2 from first: " << 0.5*pow(integrated_first/analytical_free, 2.) << "\n";
    std::cout << "0.5*I^2 from second: " << integrated_second/analytical_free << "\n";
    
    std::cout << "second_order_kft = " << second_kft << ", and integrate_epsilon_second = " << integrated_second << "\n";
    
    double theta = 0.1;
    for (int i = 0; i < 8; i++) {
        double y = i*1.0;
        std::cout << "At y=" << y << "\n";
        double epsilon_1 = epsilon_first(y, theta, mass, 3.0, rtols, atols, r_here, Mvir_over_Msun, GaussLaguerreNodes);
        std::cout << "From first order:   epsilon = " << epsilon_1 << "\n";
        
        double epsilon_2 = epsilon_second(y, theta, mass, 3.0, rtols_2, atols_2, r_here, Mvir_over_Msun, GaussLaguerreNodes, 3);
        std::cout << "From second order:  epsilon = " << sqrt(2.*epsilon_2) << "\n";
    }
    
    
    return 0;
}
