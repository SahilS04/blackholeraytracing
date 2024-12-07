#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include "boyer_lindquist_metric.h"
#include "rk45_dormand_prince.h"
#include <iomanip>
#include <omp.h>

using namespace std;
double a = 0.99;
double M = 1.0;


double intensityCalc(vector<vector<double>> result, int s, boyer_lindquist_metric metric){
    double intensity = 0;
    double r_in = 6*M;
    double r_out = 18*M;
    double r_h = M + sqrt(M*M - a*a);

    if ((r_in < abs(result[s][0]) && r_out > abs(result[s][0])) && abs(result[s][1]-(M_PI/2)) < 0.01){
        double omega = 1/(a+ (pow(result[s][0], 1.5)/sqrt(M)));
        double u_t = sqrt(metric.gamma11*pow(result[s][3], 2) + metric.gamma22*pow(result[s][4], 2) + metric.gamma33*pow(result[s][5], 2))/metric.alpha;
        double u_sub_t = -1*pow(metric.alpha, 2)*u_t - result[s][5]*metric.beta3; 
        double one_plus_z = (1 + omega*(-result[s][5]/u_sub_t))/sqrt(-metric.g_00 - omega*omega*metric.g_33 - 2*omega*metric.g_03);
        intensity = 1/pow(one_plus_z, 3);
    }

    return intensity;
}

int main(){
    double D = 500.0;
    double view_angle = (85.0/180.0)*M_PI;
    double intensity = 0;
    boyer_lindquist_metric metric = boyer_lindquist_metric(a, M);

    int num_threads = omp_get_max_threads();
    vector<rk45_dormand_prince> solvers(num_threads, rk45_dormand_prince(6, 1e-10, 1e-10));
    vector<boyer_lindquist_metric> metrics(num_threads, boyer_lindquist_metric(a, M));

    double alpha_min = -0.0230;
    double alpha_max = 0.0263;  
    double beta_min = -0.06;  
    double beta_max = 0.0623;   
    int alpha_steps = 1280;
    int beta_steps = 720;

    double alpha_step = (alpha_max - alpha_min)/alpha_steps;
    double beta_step = (beta_max - beta_min)/beta_steps; 

    double r_in = 6*M;
    double r_out = 18*M;
    double r_h = M + sqrt(M*M - a*a);
    double t0 = 0;
    double t_end = 7000;
    int size = 10000;

    vector<double> t_dense(size);
    for (int i = 0; i < size; i++){
        t_dense[i] = (t0 + i*(t_end - t0)/size);
        }

    ofstream file;
    file.open("image.csv");
    #pragma omp parallel
    {
    stringstream local_output;
    #pragma omp for
    for (int i = 0; i <= alpha_steps; i++) {
        double alpha = alpha_min + i * alpha_step;  // Use i to calculate alpha
        int thread_id = omp_get_thread_num();  // Get thread ID once per outer loop iteration
        boyer_lindquist_metric &metric = metrics[thread_id];
        rk45_dormand_prince &photon = solvers[thread_id];
        for (int j = 0; j <= beta_steps; j++) {
            intensity = 0;
            double beta = beta_min + j * beta_step;  // Use j to calculate beta
            double x_sc = D*beta;
            double y_sc = D*alpha;

            double r0 = sqrt(D*D + x_sc*x_sc + y_sc*y_sc);
            double th0 = view_angle - alpha;
            double phi0 = beta;
            metric.compute_metric(r0, th0);
            double u_r0 = -sqrt(metric.g_11)*cos(beta)*cos(alpha);
            double u_th0 = -sqrt(metric.g_22)*sin(alpha);
            double u_phi0 = sqrt(metric.g_33)*sin(beta)*cos(alpha);
            vector<double> y0 = {r0, th0, phi0, u_r0, u_th0, u_phi0};
    
            auto dydx = [&metric](double x, const vector<double> &y){
                metric.compute_metric(y[0], y[1]);
                vector<double> dydx(6);
                double ur_2 = pow(y[3], 2);
                double uth_2 = pow(y[4], 2);
                double uphi_2 = pow(y[5], 2);

                double u_t = sqrt(metric.gamma11*ur_2 + metric.gamma22*uth_2 + metric.gamma33*uphi_2)/metric.alpha;

                dydx[0] = metric.gamma11*(y[3]/u_t);
                dydx[1] = metric.gamma22*(y[4]/u_t);
                dydx[2] = (metric.gamma33*(y[5]/u_t)) - metric.beta3;
                dydx[3] = -metric.alpha*u_t*metric.d_alpha_dr + y[5]*metric.d_beta3_dr - (1/(2*u_t))*(ur_2*metric.d_gamma11_dr + uth_2*metric.d_gamma22_dr + uphi_2*metric.d_gamma33_dr);
                dydx[4] = -metric.alpha*u_t*metric.d_alpha_dth + y[5]*metric.d_beta3_dth - (1/(2*u_t))*(ur_2*metric.d_gamma11_dth + uth_2*metric.d_gamma22_dth + uphi_2*metric.d_gamma33_dth);
                dydx[5] = 0;

                return dydx;
            };

            auto stop = [t_end, r_in, r_out, r_h] (double x, vector<double> &y){
                if (y[0] <= 1.02*r_h){
                    //cout << "Hit horizon" << endl;
                    return true;
                } else if ((r_in < y[0] && r_out > y[0]) && abs(y[1]-(M_PI/2)) < 0.01){
                    //cout << "Hit disk" << endl;
                    return true;
                } else if (y[0] > 1e3){
                    return true;
                } else{
                    return x >= t_end;
                }
            };
            //cout << i << " " << j << endl;

            photon.integrate(dydx, stop, t0, y0, false, t_dense);
            int end = photon.xs.size();
             local_output << x_sc << "," << y_sc << "," << intensityCalc(photon.result, end-1, metric) << endl;
            }
        }

        #pragma omp critical
        {
            file << local_output.str();
        }
    }


    return 0;
}