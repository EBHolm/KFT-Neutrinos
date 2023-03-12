//
//  quadrature.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include "quadrature.hpp"

/*
    Computes nodes and weights for generalized Gauss-Laguerre quadrature
 
    Works up to about 90 nodes
 
    Adapted from "compute_Laguerre" of CLASS
    https://github.com/AarhusCosmology/CLASSpp/blob/master/tools/quadrature.c
 */
std::tuple<std::vector<double>, std::vector<double>> ComputeLaguerre(int N, double alpha) {
    int maxiter = 10;
    double eps = 1e-14; // tolerance of root
    
    std::vector<double> x(N);
    std::vector<double> x_weights(N);
    std::vector<double> b(N);
    std::vector<double> c(N);
    
    for(int i = 0; i < N; i++) {
        b[i] = alpha + 2.0*i + 1.0;
        c[i] = i*(alpha+i);
    }
    
    double logprod = 0.0;
    for(int i = 1; i < N; i++) {
        logprod += log(c[i]);
    }
    double logcc = lgamma(alpha + 1) + logprod;
    
    double x0 = 0.0;
    double p1 = 0.0;
    double dp2 = 0.0;
    
    /* Loop over roots: */
    for (int i = 0; i < N; i++) {
        /* Estimate root: */
        if (i == 0) {
            x0 = (1.0 + alpha)*(3.0 + 0.92*alpha)/(1.0 + 2.4*N + 1.8*alpha);
        }
        else if (i == 1) {
            x0 += (15.0 + 6.25*alpha)/(1.0 + 0.9*alpha + 2.5*N);
        }
        else {
            double r1 = (1.0+2.55*(i-1))/( 1.9*(i-1));
            double r2 = 1.26*(i-1)*alpha/(1.0+3.5*(i-1));
            double ratio = (r1+r2)/(1.0+0.3*alpha);
            x0 += ratio*(x0 - x[i-2]);
        }
        /* Refine root using Newtons method: */
        for (int iter = 1; iter <= maxiter; iter++) {
            /* We need to find p2=L_N(x0), dp2=L'_N(x0) and
               p1 = L_(N-1)(x0): */
            p1 = 1.0;
            double dp1 = 0.0;
            double p2 = x0 - alpha - 1.0;
            dp2 = 1.0;
            for (int j = 1; j < N; j++) {
                double p0 = p1;
                double dp0 = dp1;
                p1  = p2;
                dp1 = dp2;
                p2  = (x0-b[j])*p1 - c[j]*p0;
                dp2 = (x0-b[j])*dp1 + p1 - c[j]*dp0;
            }
            /* New guess at root: */
            double d = p2/dp2;
            x0 -= d;
            if (fabs(d) <= eps*(fabs(x0) + 1.0)) {
                break;
            }
        }
        /* Okay, write root and weight: */
        x[i] = x0;

        /* "totalweight" input from CLASS always true in our use-case */
        // x_weights[i] = exp(logcc-log(dp2*p1)); // These weights match those from numpy.polynomial.laguerre.laggauss
        x_weights[i] = exp(x0+logcc-log(dp2*p1)); // But it is smarter to incorporate the exp(x0) factor here
    };
    
    /*
    std::cout << "Created Gauss-Laguerre grid with the following nodes:\n";
    for (int i = 0; i < N; i++) {
        std::cout << "x=" << x[i] << ", weight=" << x_weights[i] << "\n";
    } */
    
    
    return std::make_tuple(x, x_weights);
}
