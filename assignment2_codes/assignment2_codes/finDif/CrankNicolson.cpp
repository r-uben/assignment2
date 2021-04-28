//
//  CrankNicolson.cpp
//  assignment2_codes
//
//  Created  on 28/4/21.
//

/// Header of the explicit difference method
#include "CrankNicolson.h"

/// Header with some important and usefule functions
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

/// First theclare the parameters for the problem, local grid and variables and the vectors required
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
    m_S0    = S0;
    m_Smax  = Smax;
    m_J     = J;
    m_I     = I;
    m_dS    = m_Smax / J;
    m_dt    = m_T / I;
    m_jStar = m_S0/m_dS;
}
void
CN::convertibleBond(bool lu, double tol, double omega)
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
        double t = i*m_dt;
        /// BOUNDARY CONDITIONS
        // Declare vectors for matrix equations
        vector<double> a = {0.}, b = {1.}, c = {0.}, d = {m_F * AUX::discountFactor(m_r, m_T-t)};
        // LU method
        vector<double> beta = {b[0]}, D = {d[0]};
        
        /// SET UP MATRIX EQUATIONS
        // PRINT_4DATA_LINE("a", "b", "c", "d")
        // PRINT_4DATA_LINE(a[0], b[0],c[0], d[0])
        for(int j=1;j< m_J;j++)
        {
            a.push_back( aFunc(j) );
            b.push_back( bFunc(j) );
            c.push_back( cFunc(j) );
            d.push_back( dFunc(j, vOld) );
            // PRINT_4DATA_LINE(aFunc(j), bFunc(j), cFunc(j), dFunc(j, vOld));
            
            // LU method
            if (lu == true)
            {
                beta.push_back( betaFunc(j, beta[j-1]) );
                D.push_back( DFunc(j, beta[j-1], d[j], D[j-1]) );
            }
        }
        // Boundary conditions at S = Smax
        a.push_back(0.);
        b.push_back(1.);
        c.push_back(0.);
        d.push_back(m_R * m_Smax);
        
        // LU method
        if (lu == true)
        {
            beta.push_back( betaFunc(m_J, beta[m_J-1]) );
            D.push_back( DFunc(m_J, beta[m_J-1], d[m_J], D[m_J-1]) );
            // Solve matrix equations with LU method
            vNew[0] = m_X;
            vNew[m_J] = D[m_J] / beta[m_J];
            for (int j=m_J-1; j >=0; j--)
                vNew[j] = prevV(j, beta, D, vNew);
        }
        else
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
              cout << " ( " << vNew[0] << " , ";
            }
            // 0 < j < jMax
            for(int j=1; j<m_J; j++)
            {
              double y = (d[j] - a[j] * vNew[j-1] - c[j]*vNew[j+1]) / b[j];
              vNew[j] = vNew[j] + omega * (y-vNew[j]);
              cout << vNew[j] << " , ";
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
                PRINT_DATA_LINE("NOT CONVERGED")
        }
        // Set old=new
        vOld = vNew;
    }
    // Finish looping through time levels
    // output the estimated option price
    double optionValue = approxPrice(vNew, S);
    // output the estimated option price
    PRINT_2DATA_LINE("Value of the option", optionValue);
}


// PDE COEFFICIENTS
double
CN::aFunc(int j)
{
    return 0.25 * (m_sigma * m_sigma * j * j - m_r * j);
}
double
CN::bFunc(int j)
{
    return - 0.5 * m_sigma * m_sigma * j * j - 0.5 * m_r - 1. / m_dt;
}
double
CN::cFunc(int j)
{
    return 0.25 * (m_sigma * m_sigma * j * j + m_r * j);
}
double
CN::dFunc(int j, vector<double> &v)
{
    double a = aFunc(j);
    double b = bFunc(j);
    double c = cFunc(j);
    double d =  - ( a * v[j-1] + (b + 2 / m_dt) * v[j] + c * v[j+1] );
    return d;
}
// LU COEFFICNETS
double
CN::betaFunc(int j, double prevBeta)
{
    return bFunc(j) - (aFunc(j) * cFunc(j-1)) / prevBeta;
}
double
CN::DFunc(int j, double prevBeta, double d, double prevD)
{
    return d - (aFunc(j) * prevD) / prevBeta;
}
double
CN::prevV(int j, vector<double> &beta, vector<double> &D, vector<double> &V)
{
    return 1 / beta[j] * (D[j] - cFunc(j) * V[j+1]);
}
// USEFUL FUNCTIONS

//double
//CN::interiorPoints(vector<double> &v, int j)
//{
//    double numerator    = A(j) * v[j+1] + B(j) * v[j] + C(j) * v[j-1];
//    double denominator  = 1. / (1. + m_r * m_dt);
//    return numerator * denominator;
//}

double
CN::approxPrice(vector<double> &v, vector<double> &s)
{
    vector<double> p1 = {s[m_jStar], v[m_jStar]};
    vector<double> p2 = {s[m_jStar+1], v[m_jStar+1]};
    return AUX::lagInterp(m_S0, p1, p2);
}
