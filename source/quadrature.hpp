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

/*
    Adaptive Gauss-Kronrod quadrature
    15-point Gauss-Kronrod and 7-point Gauss quadratures
*/

template <class FunctionArguments>
std::tuple<double, double> GaussKronrod_quad(std::function<double(double, FunctionArguments)> func, double a, double b, FunctionArguments f_args) {
    std::vector<double> z_k = {-0.991455371120813, -0.949107912342759, -0.864864423359769, -0.741531185599394,
                               -0.586087235467691, -0.405845151377397, -0.207784955007898,  0.0,
                                0.207784955007898,  0.405845151377397,  0.586087235467691,  0.741531185599394,
                                0.864864423359769,  0.949107912342759,  0.991455371120813};
    std::vector<double> w_k = { 0.022935322010529,  0.063092092629979,  0.104790010322250,  0.140653259715525,
                                0.169004726639267,  0.190350578064785,  0.204432940075298,  0.209482141084728,
                                0.204432940075298,  0.190350578064785,  0.169004726639267,  0.140653259715525,
                                0.104790010322250,  0.063092092629979,  0.022935322010529};
    std::vector<double> w_g = { 0.129484966168870,  0.279705391489277,  0.381830050505119,  0.417959183673469,
                                0.381830050505119,  0.279705391489277,  0.129484966168870};
    
    std::vector<double> y(15);
    for (int i = 0; i < 15; i++) {
        double t = 0.5*(a*(1 - z_k[i]) + b*(1 + z_k[i]));
        y[i] = func(t, f_args);
    }
    double I_k = 0.0, I_g = 0.0;
    for (int i_k = 0; i_k < 15; i_k++) {
        /* Kronrod quadrature */
        w_k[i_k] *= 0.5*(b - a);
        I_k += y[i_k]*w_k[i_k];
    }
    for (int i_g = 0; i_g < 7; i_g++) {
        /* Gauss quadrature */
        w_g[i_g] *= 0.5*(b - a);
        I_g += y[2*i_g + 1]*w_g[i_g];
    }
    double err = pow(200.0*fabs(I_k - I_g), 1.5);
    return std::make_tuple(I_k, err);
}

template <class FunctionArguments>
std::tuple<double, double> GaussKronrod(std::function<double(double, FunctionArguments)> func, double a, double b, double rtol, double atol, FunctionArguments f_args) {
    auto [I, err] = GaussKronrod_quad(func, a, b, f_args);
    
    if (((fabs(I) > 0 ) && (err/I < rtol)) || (err < atol)) {
        return std::make_tuple(I, err);
    }
    else {
        double m = 0.5*(a + b);
        auto [I_left, err_left]   = GaussKronrod(func, a, m, 1.5*rtol, atol, f_args);
        auto [I_right, err_right] = GaussKronrod(func, m, b, 1.5*rtol, atol, f_args);
        return std::make_tuple(I_left + I_right, sqrt(pow(err_left, 2.) + pow(err_right, 2.)));
    }
}


std::tuple<std::vector<double>, std::vector<double>> ComputeLaguerre(int N, double alpha = 0.0);


template <class FunctionArguments>
double GaussLaguerre(std::function<double(double, FunctionArguments)> func, int N, FunctionArguments args) {
    auto [x, x_weights] = ComputeLaguerre(N);
    double val = 0.0;
    for (auto i = 0; i < x.size(); i++) {
        val += func(x[i], args)*x_weights[i];
    }
    return val;
}



#endif /* quadrature_hpp */
