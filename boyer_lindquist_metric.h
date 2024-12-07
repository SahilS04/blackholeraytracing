#include <iostream>
#include <cmath>
#include <omp.h>
#include <vector>
#pragma once

using namespace std;

class boyer_lindquist_metric {
public:
    boyer_lindquist_metric(double a0, double M0) {
        a = a0;
        M = M0;
    }

void cout_metric(){
    cout << "rho: " << rho << endl;
    cout << "delta: " << delta << endl;
    cout << "sigma: " << sigma << endl;
    cout << "rho_prime_r: " << rho_prime_r << endl;
    cout << "delta_prime_r: " << delta_prime_r << endl;
    cout << "sigma_prime_r: " << sigma_prime_r << endl;
    cout << "rho_prime_th: " << rho_prime_th << endl;
    cout << "sigma_prime_th: " << sigma_prime_th << endl;
    cout << "alpha: " << alpha << endl;
    cout << "beta3: " << beta3 << endl;
    cout << "gamma11: " << gamma11 << endl;
    cout << "gamma22: " << gamma22 << endl;
    cout << "gamma33: " << gamma33 << endl;
    cout << "g_00: " << g_00 << endl;
    cout << "g_03: " << g_03 << endl;
    cout << "g_11: " << g_11 << endl;
    cout << "g_22: " << g_22 << endl;
    cout << "g_33: " << g_33 << endl;
    cout << "d_alpha_dr: " << d_alpha_dr << endl;
    cout << "d_beta3_dr: " << d_beta3_dr << endl;
    cout << "d_gamma11_dr: " << d_gamma11_dr << endl;
    cout << "d_gamma22_dr: " << d_gamma22_dr << endl;
    cout << "d_gamma33_dr: " << d_gamma33_dr << endl;
    cout << "d_alpha_dth: " << d_alpha_dth << endl;
    cout << "d_beta3_dth: " << d_beta3_dth << endl;
    cout << "d_gamma11_dth: " << d_gamma11_dth << endl;
    cout << "d_gamma22_dth: " << d_gamma22_dth << endl;
    cout << "d_gamma33_dth: " << d_gamma33_dth << endl;
}


void compute_metric(double r, double th) {
    double episilon = 1.0e-12;

    double sint = sin(th);
    double cost = cos(th);
    double sint2 = sint*sint;
    double cost2 = cost*cost;
    double r2 = r*r;
    double a2 = a*a;

    rho = r2 + a2 * cost2;
    delta = r2 - 2.0*M*r + a2;
    sigma = pow(r2+a2, 2) - (a2 * delta * sint2);

    rho_prime_r = 2.0*r;
    delta_prime_r = 2.0*r - 2.0*M;
    sigma_prime_r = 4.0*r*(r2+a2) - a2*sint2*delta_prime_r;

    rho_prime_th = -2.0*a2*cost*sint;
    sigma_prime_th = 2.0*a2*delta*cost*sint;

    alpha = sqrt((rho*delta)/sigma);
    beta3 = (-2.0*M*r*a)/sigma;
    gamma11 = delta/rho;
    gamma22 = 1.0/rho;
    gamma33 = rho/(sigma*sint2);
    g_00 = ((2*M*r)/rho) - 1.0;
    g_03 = (-2*M*a*r*sint2)/rho;
    g_11 = rho/delta; // divide by zero error
    g_22 = rho;
    g_33 = (sigma*sint2)/rho;

    double rho2 = pow(rho, 2);
    double sigma2 = pow(sigma, 2);
    
    d_alpha_dr = (((rho_prime_r*delta)+(rho*delta_prime_r))*sigma - sigma_prime_r*rho*delta)/(2*sigma2*alpha);
    d_alpha_dth = ((rho_prime_th*delta*sigma) - sigma_prime_th*rho*delta)/(2*sigma2*alpha);

    d_beta3_dr = (-2*M*a*sigma + 2*M*a*r*sigma_prime_r)/sigma2;
    d_beta3_dth = (2*M*a*r*sigma_prime_th)/sigma2;

    d_gamma11_dr = (delta_prime_r*rho - rho_prime_r*delta)/rho2;
    d_gamma11_dth = (-rho_prime_th*delta)/rho2;

    d_gamma22_dr = (-rho_prime_r)/rho2;
    d_gamma22_dth = (-rho_prime_th)/rho2;

    d_gamma33_dr = (rho_prime_r*sigma - sigma_prime_r*rho)/(sigma2 * sint2);
    d_gamma33_dth = (rho_prime_th*sigma*sint - (sigma_prime_th*sint + 2*sigma*cost)*rho)/(sigma2*pow(sint, 3));
}


double a, M;
double alpha, beta3;
double rho, delta, sigma;
double delta_prime_r, sigma_prime_r, sigma_prime_th;
double rho_prime_th, rho_prime_r;
double gamma11, gamma22, gamma33; // components of upper gamma^ij
double g_00, g_03, g_11, g_22, g_33; // components of lower g_\mu\nu
double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr;
double d_alpha_dth, d_beta3_dth, d_gamma11_dth, d_gamma22_dth, d_gamma33_dth;

};