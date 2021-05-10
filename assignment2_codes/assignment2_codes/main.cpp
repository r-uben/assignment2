#include <cmath>
#include <vector>
#include <algorithm>

#include <iostream>
#include <fstream>
// For this task
#include "CrankNicolson.h"
#include "ConvertibleBonds.h"
#include "Q1.h"
#include "Q2.h"
// General
#include "GeneralFunctions.h"
// To calculate the time
#include <chrono>
#define  CHRONO   std::chrono
#define  SET_TIME CHRONO::system_clock::now()
#define  DURATION CHRONO::duration
#define  MILLI    std::milli

#define NEXT_CASE cout << "#######" << endl;

using namespace std;
/* Template code for the Crank Nicolson Finite Difference
 */
int main(){
    cout.precision(15);
    // PARAMETERS
    double T = 3., F = 35., R = 2, r = 0.0117, kappa = 0.0833333333333, mu = 0.0053, X = 17.38, C = 0.205, alpha = 0.02, beta = 0.808, sigma = 0.66;
    double Smax = 250;
    double I = 500;
    double J = 500;
    ofstream output;
    Q1 q1(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 1000, I, J);
    Q2 q2(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 1000, I, J);
    /*
        TASK 1
    */
    //    CONV_BONDS convBond(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, Smax, J, I);
    //    cout << (F - convBond.B(0)) / convBond.A(0) << endl;
    
    ///
    /// FIXING BETA AND MOVING SIGMA
    ///
        // vector <double> sigmas = {0, 0.381, 0.66, 0.802, 0.962, 1};
        // q1.variousSigmas(sigmas, THOMAS);
    
    ///
    /// FIXING SIGMA AND MOVING BETA
    ///
        // vector <double> betas = {0.892};
        // q1.variousBetas(betas, THOMAS);
    
    ///
    /// FIXING SIGMA AND AND BETA AND MOVING KAPPA
    ///
        // vector <double> kappas = {0, 0.02, 0.0833, 0.2, 0.5, 1};
        // q1.variousKappas(kappas, THOMAS);
    
    ///
    /// INCREASING IMAX AND JMAX (SQUARE) AND FIXING SMAX
    ///
        //  q1.variousV_fixedS0(17.38, 2, 2, 20000, 3, 21.);
        //  vector<double> v = {42.0495722630, 42.0495647830, 42.0495620065};
        //  q1.V_fixedS0(17.38, 16, 8., 952, 952);
        //  q1.V_fixedS0(17.38, 16, 8., 9216, 9216);
        //  AUX::extrap(v[1], v[2], 2);

    /*
        TASK 2
    */
    //    q1.increasingS(500, 500, 300, 1.05);
    //    q2.increasingS(100, 100, 200, 1.05);
    //    q2.increasingS(200, 200, 200, 1.05);
    //    vector<double> rs = {0.00585, 0.0117, 0.01755};
    //    q2.variousInterestRates(rs, 300, 300, 200);
          q2.variousV_fixedS0(17.38, 0.96151, 2, 2, 80000, 16, 8.);
    //    q2.V_fixedS0(17.38, 0.96151, 16, 8., 128, 64);
    //    q2.V_fixedS0(17.38, 0.96151, 16, 8., 256, 128);
    //    AUX::extrap(43.18821969729999, 43.2219325264, 2);
    return 0;
}
