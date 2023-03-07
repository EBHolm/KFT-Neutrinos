//
//  quadrature.cpp
//  KFT-Neutrinos
//
//  Created by Emil Brinch Holm on 01/03/2023.
//

#include "quadrature.hpp"


/*
    Adaptive Gauss-Kronrod quadrature
    15-point Gauss-Kronrod and 7-point Gauss quadratures
*/

//template <class FunctionArguments>
std::tuple<double, double> GaussKronrod_quad(std::function<double(double, FirstOrderArguments)> func, double a, double b, FirstOrderArguments f_args) {
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

//template <class FunctionArguments>
std::tuple<double, double> GaussKronrod(std::function<double(double, FirstOrderArguments)> func, double a, double b, double rtol, double atol, FirstOrderArguments f_args) {
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

/*
    Gauss-Laguerre quadrature
 
*/

//template <class FunctionArguments>
double GaussLaguerre(std::function<double(double, FirstOrderArguments)> func, int N, FirstOrderArguments args) {
    auto [x, x_weights] = ComputeLaguerre(N);
    double val = 0.0;
    for (auto i = 0; i < x.size(); i++) {
        val += func(x[i], args)*x_weights[i];
    }
    return val;
}

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


// ------------ Warning: Code duplication, second order versions of the same ---------------------------------

std::tuple<double, double> GaussKronrod_quad(std::function<double(double, SecondOrderArguments)> func, double a, double b, SecondOrderArguments f_args) {
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
std::tuple<double, double> GaussKronrod(std::function<double(double, SecondOrderArguments)> func, double a, double b, double rtol, double atol, SecondOrderArguments f_args) {
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
double GaussLaguerre(std::function<double(double, SecondOrderArguments)> func, int N, SecondOrderArguments args) {
    auto [x, x_weights] = ComputeLaguerre(N);
    double val = 0.0;
    for (auto i = 0; i < x.size(); i++) {
        val += func(x[i], args)*x_weights[i];
    }
    return val;
}
