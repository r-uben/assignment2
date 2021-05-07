
#include <cmath>
#include <vector>
#include <algorithm>

#include "CrankNicolson.h"
#include "ConvertibleBonds.h"
#include "Q1.h"
#include "Q2.h"


#include <iostream>
#include <fstream>

#define NEXT_CASE cout << "#######" << endl;

using namespace std;
/* Template code for the Crank Nicolson Finite Difference
 */
int main(){
    cout.precision(15);
    
    /*
        TASK 1
    */
    
    // PARAMETERS
    double T = 3., F = 35., R = 2, r = 0.0117, kappa = 0.0833333333333, mu = 0.0053, X = 17.38, C = 0.205, alpha = 0.02, beta = 0.808, sigma = 0.66;
    double Smax = 250;
    double I = 500;
    double J = 500;
    //    CONV_BONDS convBond(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, Smax, J, I);
    //    cout << (F - convBond.B(0)) / convBond.A(0) << endl;
    
    ofstream output;
    //Q1 q1(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 8*X, I, J);
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
    //    q1.fixedS0(17.38, 2, 4, 20000, 2, 4);
    //    q1.fixedS0(17.38, 2, 4, 20000, 3, 4);
    //    q1.fixedS0(17.38, 2, 4, 20000, 4, 4);
    //    q1.fixedS0(17.38, 2, 4, 20000, 8, 4);
    //    q1.fixedS0(17.38, 2, 4, 20000, 16, 4);
    //
    //    q1.fixedS0(17.38, 2, 6, 20000, 2, 6);
    //    q1.fixedS0(17.38, 2, 6, 20000, 3, 6);
    //    q1.fixedS0(17.38, 2, 6, 20000, 4, 6);
    //    q1.fixedS0(17.38, 2, 6, 20000, 8, 6);
    //    q1.fixedS0(17.38, 2, 6, 20000, 16, 6);
    //
    //    q1.fixedS0(17.38, 2, 8, 20000, 2, 8);
    //    q1.fixedS0(17.38, 2, 8, 20000, 3, 8);
    //    q1.fixedS0(17.38, 2, 8, 20000, 4, 8);
    //    q1.fixedS0(17.38, 2, 8, 20000, 8, 8);
    //    q1.fixedS0(17.38, 2, 8, 20000, 16, 8);

    
    /*
        TASK 2
    */
    
    Q1 q1(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 300*3, I, J);
    Q2 q2(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 300*3, I, J);
    q2.increasingS(500, 500, 300, 1.05);
    //q1.increasingS(500, 500, 300, 1.05);
    vector<double> rs = {0.00585, 0.0117, 0.01755};
    // q2.variousInterestRates(rs, 300, 300, 200);
    // q1.fixedS0(17.38, 2, 4, 20000, 16, 4);
    //q2.fixedS0(17.38, 2, 4, 20000, 16, 4);

}





///
/// INCREASING JMAX
///

//    Smax = 50;
//    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_increasing_jMax_iMax500_Smax50.csv");
//    double S0 = 17.38;
//    output << "F,S0,Smax,V" << endl;
//    for (int j=50; j<=5000; j=j+50)
//    {
//        CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, S0, 50, j, 500);
//        crank.eurConvertibleBond(&output, THOMAS);
//    }
//    output.close();
//    // Smax = 100
//    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_increasing_jMax_iMax500_Smax100.csv");
//    output << "F,S0,Smax,V" << endl;
//    for (int j=50; j<=5000; j=j+50)
//    {
//        CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, S0, 100, j, 500);
//        crank.eurConvertibleBond(&output, THOMAS);
//    }
//    output.close();
//    cout << "#######" << endl;
//    // Smax = 250
//    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_increasing_jMax_iMax500_Smax250.csv");
//    output << "F,S0,Smax,V" << endl;
//    for (int j=50; j<=5000; j=j+50)
//    {
//        CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, S0, 250, j, 500);
//        crank.eurConvertibleBond(&output, THOMAS);
//    }
//    output.close();
//    cout << "#######" << endl;
//    // Smax = 500
//    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_increasing_jMax_iMax500_Smax500.csv");
//    output << "F,S0,J,V" << endl;
//    for (int j=50; j<=5000; j=j+50)
//    {
//        CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, S0, 500, j, 500);
//        crank.eurConvertibleBond(&output, THOMAS);
//    }
//    output.close();

///
/// INCREASING IMAX
///
//    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_increasing_iMax_jMax5000_Smax100.csv");
//    double S0 = 17.38;
//    output << "F,S0,I,V" << endl;
//    for (int i=5; i<=500; i=i+5)
//    {
//        CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, S0, Smax, J, i);
//        crank.eurConvertibleBond(&output, THOMAS);
//    }
