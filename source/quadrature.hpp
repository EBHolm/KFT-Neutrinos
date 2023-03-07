//
//  quadrature.hpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#ifndef quadrature_hpp
#define quadrature_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "utils.hpp"

struct FirstOrderArguments { // Could make this a struct now
    double rtols[3];
    double atols[3];
    int GaussLaguerreNodes;
    
    double mnu;
    double Tnu;
    double r_here;
    double rr_here;
    double G0_Gz;
    double weight;
    double Rs;
    double p;
    double gp;
    double r_a;
    double r;
    double phi_integrand;
};

struct SecondOrderArguments { // Could make this a struct now
    double rtols[4];
    double atols[4];
    int GaussLaguerreNodes;
    
    double mnu;
    double Tnu;
    double r_here;
    double rr_here;
    
    std::vector<double> LaguerreNodes;
    std::vector<double> LaguerreWeights;
};


std::tuple<double, double> GaussKronrod(std::function<double(double, FirstOrderArguments)> func, double a, double b, double rtol, double atol, FirstOrderArguments f_args);
std::tuple<double, double> GaussKronrod(std::function<double(double, SecondOrderArguments)> func, double a, double b, double rtol, double atol, SecondOrderArguments f_args);

std::tuple<double, double> GaussKronrod_quad(std::function<double(double, FirstOrderArguments)> func, double a, double b, FirstOrderArguments f_args);
std::tuple<double, double> GaussKronrod_quad(std::function<double(double, SecondOrderArguments)> func, double a, double b, SecondOrderArguments f_args);

double GaussLaguerre(std::function<double(double, FirstOrderArguments)> func, int N, FirstOrderArguments args);
double GaussLaguerre(std::function<double(double, SecondOrderArguments)> func, int N, SecondOrderArguments args);
std::tuple<std::vector<double>, std::vector<double>> ComputeLaguerre(int N, double alpha = 0.0);

#endif /* quadrature_hpp */
