//
//  CrankNicolson.cpp
//  assignment2_codes
//
//  Created  on 28/4/21.
//

/// Header of the explicit difference method
#include "CrankNicolson.h"

/// Header with convertible bonds functions
#include "ConvertibleBonds.h"

/// Header with some important and useful functions
#include "GeneralFunctions.h"

/// Header with payoff functions of several options
#include "PayoffFunctions.h"

/// Header for constructing the grid of the explicit difference method
#include "GridConstructor.h"

/// Header for printing data
#include "PrintMacros.h"

/// Important libraries
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

/// First declare the parameters for the problem, local grid and variables and the vectors required
/// In option pricing problems, we must for S over the semi infinite domain, therefore numerically we need to choose an appropiate Smax
/// The first five parameters are intrinsice option parameters, whilst the last three ones are those appropriate parameters for the method.

CN::CCrankNicolson(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double S0, double Smax, int J, int I){
    m_T         = T;
    m_F         = F;
    m_R         = R;
    m_r         = r;
    m_kappa     = kappa;
    m_mu        = mu;
    m_X         = X;
    m_C         = C;
    m_alpha     = alpha;
    m_beta      = beta;
    m_sigma     = sigma;
    m_kappar    = kappa + r;
    m_alphar    = alpha + r;
    m_S0        = S0;
    m_Smax      = Smax;
    m_J         = J;
    m_I         = I;
    m_dS        = m_Smax / J;
    m_dt        = m_T / I;
    m_jStar     = m_S0/m_dS;
}
void
CN::convertibleBond(ofstream *output, bool lu, double tol, double omega)
{
    vector <double> vOld(m_J+1), vNew(m_J+1);
    // Setup and initialise the stock price
    vector <double> S = GRID::setupStockPrices(m_dS, m_J);
    // Setup and initialise the final conditions on the option price
    for (int j=0; j<=m_J; j++)
    {
        vOld[j] = PAYOFF::convBond(S[j], m_R, m_F);
        vNew[j] = PAYOFF::convBond(S[j], m_R, m_F);
        // PRINT_5DATA_LINE(m_I, j, S[j], vNew[j], vOld[j])
    }
    // start looping through time levels
    for(int i=m_I-1; i>=0; i--)
    {
        // CONVERTIBLE BONDS LIBRARY
        CONV_BONDS convBonds(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma);
        double t = i*m_dt;
        /// BOUNDARY CONDITIONS
        // Declare vectors for matrix equations
        vector<double> a = {0.}, b = {1.}, c = {0.}, d = {convBonds.V_S0(t)};
        // LU method
        vector<double> beta = {b[0]}, D = {d[0]};
        
        /// SET UP MATRIX EQUATIONS
        // PRINT_4DATA_LINE("a", "b", "c", "d")
        // PRINT_4DATA_LINE(a[0], b[0],c[0], d[0])
        for(int j=1;j< m_J;j++)
        {
            a.push_back(aFunc(t, j));
            b.push_back(bFunc(t, j));
            c.push_back(cFunc(t, j));
            d.push_back(dFunc(t, j, vOld));
            
            // LU method
            if (lu == true)
            {
                beta.push_back(betaFunc(t, j, beta[j-1]));
                D.push_back(DFunc(t, j, beta[j-1], d[j], D[j-1]));
            }
        }
        // Boundary conditions at S = Smax
        a.push_back(0.);
        b.push_back(1.);
        c.push_back(0.);
        d.push_back(convBonds.V(m_Smax,t));
        
        // LU METHOD
        if (lu == true)
        {
            beta.push_back( betaFunc(t, m_J, beta[m_J-1]) );
            D.push_back( DFunc(t, m_J, beta[m_J-1], d[m_J], D[m_J-1]) );
            // Solve matrix equations with LU method
            vNew[0] = m_X;
            vNew[m_J] = D[m_J] / beta[m_J];
            for (int j=m_J-1; j >=0; j--)
                vNew[j] = prevV(t, j, beta, D, vNew);
        }
        if (lu == false)
        {
            // Solve matrix equations with SOR
            int sor, iterMax = 10000;
            for (sor = 0; sor < iterMax; sor++)
            {
            cout << sor << " == \n";
            // SOR equations in here
            // j = 0
            {
              double y = (d[0] - c[0] * vNew[1]) / b[0];
              vNew[0] = vNew[0] + omega * (y-vNew[0]);
              cout << " (" << vNew[0] << ", ";
            }
            // 0 < j < jMax
            for(int j=1; j<m_J; j++)
            {
              double y = (d[j] - a[j] * vNew[j-1] - c[j]*vNew[j+1]) / b[j];
              vNew[j] = vNew[j] + omega * (y-vNew[j]);
              cout << vNew[j] << ", ";
            }
            // j = jMax
            {
              double y = (d[m_J] - a[m_J] * vNew[m_J-1]) / b[m_J];
              vNew[m_J] = vNew[m_J] + omega * (y-vNew[m_J]);
              cout << vNew[m_J] << " )\n";
            }
            // calculate residual
            double error=0.;
            error += fabs(d[0] - b[0] * vNew[0] - c[0] * vNew[1]);

            for(int j=1; j<m_J;j++)
              error += fabs(d[j] - a[j]*vNew[j-1] - b[j]*vNew[j] - c[j]*vNew[j+1]);

            error += fabs(d[m_J] - a[m_J]*vNew[m_J-1] - b[m_J]*vNew[m_J]);
            // check for convergence and exit loop if converged
            if(error<tol)
              break;
            }
            if(sor==iterMax)
                PRINT_DATA_LINE("NOT CONVERGED");
        }
        // Set old=new
        vOld = vNew;
    }
    // Finish looping through time levels
    // output the estimated option price
    double optionValue = approxPrice(vNew, S);
    // output the estimated option price
    //ofstream output;
    // OpenCSVFile(&output, "eurConvBondValues", NOT_OVER_WRITE);
    DATA_LINE(output, m_F, m_S0, optionValue);
    PRINT_DATA_LINE(m_S0, optionValue);
}


// PDE COEFFICIENTS
double
CN::aFunc(double t, int j)
{
    double first_term   = -0.25 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2*(m_beta-1));
    double second_term  =  0.25 * m_kappa * ( theta(t) / m_dS - j);
    return first_term + second_term;
}

double
CN::bFunc(double t, int j)
{
    double long_term = 0.5 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2.*(m_beta-1.));
    return 1. / m_dt + 0.5 * m_r  + long_term;
}

double
CN::cFunc(double t, int j)
{
    double first_term   = -0.25 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2.*(m_beta-1.));
    double second_term  = -0.25 * m_kappa * ( theta(t) / m_dS - j);
    return first_term + second_term;
}

double
CN::dFunc(double t, int j, vector<double> &v)
{
    double a = aFunc(t, j);
    double b = bFunc(t, j);
    double c = cFunc(t, j);
    double d =  - ( a * v[j-1] + (b - 2 / m_dt) * v[j] + c * v[j+1] ) + m_C * exp(-m_alpha * t);
    return d;
}


// LU COEFFICIENTS
double
CN::betaFunc(double t, int j, double prevBeta)
{
    return bFunc(t, j) - (aFunc(t, j) * cFunc(t, j-1)) / prevBeta;
}
double
CN::DFunc(double t, int j, double prevBeta, double d, double prevD)
{
    return d - (aFunc(t, j) * prevD) / prevBeta;
}
double
CN::prevV(double t, int j, vector<double> &beta, vector<double> &D, vector<double> &V)
{
    return 1 / beta[j] * (D[j] - cFunc(t, j) * V[j+1]);
}
// USEFUL FUNCTIONS
double
CN::theta(double t)
{
    return (1 + m_mu) * m_X * exp(m_mu * t);
}

double
CN::approxPrice(vector<double> &v, vector<double> &s)
{
    vector<double> p1 = {s[m_jStar], v[m_jStar]};
    vector<double> p2 = {s[m_jStar+1], v[m_jStar+1]};
    return AUX::lagInterp(m_S0, p1, p2);
}
