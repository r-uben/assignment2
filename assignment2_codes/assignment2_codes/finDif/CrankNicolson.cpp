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

// To calculate the time
#include <chrono>
#define  CHRONO   std::chrono
#define  SET_TIME CHRONO::system_clock::now()
#define  DURATION CHRONO::duration
#define  MILLI    std::milli

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
    m_S0    = S0;
    m_Smax  = Smax;
    m_J     = J;
    m_I     = I;
    m_dS    = m_Smax / J;
    m_dt    = m_T / I;
    m_jStar = m_S0/m_dS;
}
void
CN::eurConvertibleBond(ofstream *output, int method, int degree, double tol, double omega)
{
    // CONVERTIBLE BONDS LIBRARY
    CONV_BONDS convBonds(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, m_Smax, m_J, m_I);
    //
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
    for(long i=m_I-1; i>=0; i--)
    {
        double t        = i * m_dt;
        double approx_t = (i + 0.5) * m_dt;
        /// BOUNDARY CONDITIONS
        // Declare vectors for matrix equations
        vector<double> a = {0.}, b = {1.}, c = {0.}, d = {convBonds.V_S0(approx_t)};
        // LU method
        vector<double> beta = {b[0]}, D = {d[0]};
        
        /// SET UP MATRIX EQUATIONS
        // PRINT_4DATA_LINE("a", "b", "c", "d")
        // PRINT_4DATA_LINE(a[0], b[0],c[0], d[0])
        for(int j=1;j< m_J;j++)
        {
            a.push_back(convBonds.aFunc(i, j));
            b.push_back(convBonds.bFunc(i, j));
            c.push_back(convBonds.cFunc(i, j));
            d.push_back(convBonds.dFunc(i, j, vOld));
            
            // LU method
            if (method == LU)
            {
                beta.push_back(convBonds.betaFunc(t, j, beta[j-1]));
                D.push_back(convBonds.DFunc(t, j, beta[j-1], d[j], D[j-1]));
            }
        }
        // Boundary conditions at S = Smax
        a.push_back(0.);
        b.push_back(1.);
        c.push_back(0.);
        d.push_back(convBonds.V_Smax(m_Smax,approx_t));
        
        // EUROPEAN OPTIONS
        // SOLVE WITH LU METHOD
        if (method == LU)
        {
            beta.push_back(convBonds.betaFunc(i, m_J, beta[m_J-1]));
            D.push_back(convBonds.DFunc(i, m_J, beta[m_J-1], d[m_J], D[m_J-1]));
            // Solve matrix equations with LU method
            vNew[0] = m_X;
            vNew[m_J] = D[m_J] / beta[m_J];
            for (long j=m_J-1; j >=0; j--)
                vNew[j] = convBonds.prevV(t, j, beta, D, vNew);
        }
        // SOLVE WITH SOR METHOD
        if (method == SOR)
        {
            // Solve matrix equations with SOR
            long y, sor, iterMax = 10000;
            for (sor = 0; sor < iterMax; sor++)
            {
                // SOR equations in here
                // j = 0
                {
                  y = (d[0] - c[0] * vNew[1]) / b[0];
                  vNew[0] = vNew[0] + omega * (y-vNew[0]);
                }
                // 0 < j < jMax
                for(long j=1; j<m_J; j++)
                {
                  y = (d[j] - a[j] * vNew[j-1] - c[j]*vNew[j+1]) / b[j];
                  vNew[j] = vNew[j] + omega * (y-vNew[j]);
                }
                // j = jMax
                {
                  y = (d[m_J] - a[m_J] * vNew[m_J-1]) / b[m_J];
                  vNew[m_J] = vNew[m_J] + omega * (y-vNew[m_J]);
                }
                // Calculate residual
                double error=0.;
                error += fabs(d[0] - b[0] * vNew[0] - c[0] * vNew[1]);
                for(long j=1; j<m_J;j++)
                  error += fabs(d[j] - a[j]*vNew[j-1] - b[j]*vNew[j] - c[j]*vNew[j+1]);
                error += fabs(d[m_J] - a[m_J]*vNew[m_J-1] - b[m_J]*vNew[m_J]);
                // Check for convergence and exit loop if converged
                if(error < tol)
                  break;
            }
            if(sor==iterMax)
                PRINT_DATA_LINE("NOT CONVERGED");
        }
            // SOLVE WITH THOMAS SOLVER ALGORITHM
            if (method == THOMAS)
                vNew = thomasSolve(a, b, c, d);
        // Set old=new
        vOld = vNew;
    }
    // Finish looping through time levels
    // output the estimated option price
    double optionValue = GRID::lagrangeInterpolation(vNew, S, m_S0, degree);//approxPrice(vNew, S);
    // output the estimated option price
    double asympV = m_S0*convBonds.A(0) + convBonds.B(0);
    DATA_LINE(output, m_F, m_I, m_J, m_S0, optionValue, asympV);
    PRINT_DATA_LINE(m_F, m_I, m_J, m_S0, optionValue, asympV);

}

void
CN::amConvertibleBond(ofstream *output, int method, int degree, double tol, double omega)
{
    //    cout << "Price for the put option: ";
    //    cin >> m_P;
    //    cout << "Until what time the option is early exercisable?: ";
    //    cin >> m_t0;
    m_P = 38.;
    m_t0 = 0.96151;
    m_iterMax = 10000;
    m_rho = 1e8;
    //    if (method == PENALTY)
    //    {
    //        cout << "Number of iterations of the Penalty Method: ";
    //        cin >> m_iterMax;
    //        cout << "Penalty parameter: ";
    //        cin >> m_rho;
    //    }
    
    // CONVERTIBLE BONDS LIBRARY
    CONV_BONDS convBonds(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, m_Smax, m_J, m_I);
    //
    vector <double> vOld(m_J+1), vNew(m_J+1);
    // Setup and initialise the stock price
    vector <double> S = GRID::setupStockPrices(m_dS, m_J);
    // Setup and initialise the final conditions on the option price
    for (long j=0; j<=m_J; j++)
    {
        vOld[j] = PAYOFF::convBond(S[j], m_R, m_F);
        vNew[j] = PAYOFF::convBond(S[j], m_R, m_F);
        // PRINT_5DATA_LINE(m_I, j, S[j], vNew[j], vOld[j])
    }
    // start looping through time levels
    for(long i=m_I-1; i>=0; i--)
    {
        double t        = i*m_dt;
        double approx_t = (i + 0.5) * m_dt;
        /// BOUNDARY CONDITIONS
        // Declare vectors for matrix equations
        vector<double> a = {0.}, b = {1.}, c = {0.}, d;
        if (t <= m_t0)
            d = {max(m_P,convBonds.V_S0(approx_t))};
        else
            d = {convBonds.V_S0(approx_t)};
        /// SET UP MATRIX EQUATIONS
        for(int j=1;j< m_J;j++)
        {
            a.push_back(convBonds.aFunc(i, j));
            b.push_back(convBonds.bFunc(i, j));
            c.push_back(convBonds.cFunc(i, j));
            d.push_back(convBonds.dFunc(i, j, vOld));
        }
        // Boundary conditions at S = Smax
        a.push_back(0.);
        b.push_back(1.);
        c.push_back(0.);
        d.push_back(m_R*m_Smax);
        // Temporal restriction
        if (method == PENALTY)
        {
            int penaltyIt;
            for (int penaltyIt=0; penaltyIt < m_iterMax; penaltyIt++)
            {
                // Create new vectors containing a copy of the FD approx
                vector<double> aHat(a), bHat(b), cHat(c), dHat(d);
                // Apply penalty here to finite difference scheme
                for (int j=1; j<m_J; j++)
                {
                    // If current value suggesta apply penalty, adjust matrix equations
                    if( vNew[j] < m_R*S[j] )
                    {
                       bHat[j] = b[j] + m_rho;
                       dHat[j] = d[j] + m_rho*(m_R*S[j]);
                    }
                    if( t <= m_t0)
                    {
                       if ( vNew[j] < m_P)
                       {
                           bHat[j] = bHat[j] + m_rho;
                           dHat[j] = dHat[j] + m_rho*m_P;
                       }
                   }
                    
                }
                // Solve with thomas Method
                vector <double> y = thomasSolve(aHat, bHat, cHat, dHat);
                // y now contains next guess at solution
                // Check for differences between vNew and y
                double error = 0.;
                for (int j=0; j<= m_J; j++)
                {
                    error += (vNew[j] - y[j])*(vNew[j] - y[j]);
                }
                // Update Value of vNew
                vNew = y;
                // make an exit condition when solution is converged
                if (error < tol * tol)
                {
                    // START_LINE "Solved after " << penaltyIt << " iterations" END_LINE;
                    break;
                }
            }
            if(penaltyIt >= m_iterMax)
            {
                PRINT_DATA_LINE("Error NOT converging within required iterations");
                throw;
            }
            vOld=vNew;
        }  // Ending penalty method
    } // Ending temporal loop
    
    // Finish looping through time levels
    // output the estimated option price
    double optionValue = GRID::lagrangeInterpolation(vNew, S, m_S0, degree);//approxPrice(vNew, S);
    // output the estimated option price
    double asympV = m_S0*convBonds.A(0) + convBonds.B(0);
    DATA_LINE(output, m_F, m_I, m_J, m_S0, optionValue, asympV);
    PRINT_DATA_LINE(m_F, m_I, m_J, m_S0, optionValue, asympV);
}


// USEFUL FUNCTIONS
double
CN::approxPrice(vector<double> &v, vector<double> &s)
{
    vector<double> p1 = {s[m_jStar], v[m_jStar]};
    vector<double> p2 = {s[m_jStar+1], v[m_jStar+1]};
    return AUX::lagInterp(m_S0, p1, p2);
}

vector<double>
CN::thomasSolve(const vector<double> &a,const vector<double> &b_,const vector<double> &c, vector<double> &d)
{
    auto n=a.size();
    std::vector<double> b(n),temp(n);
    // initial first value of b
    b[0]=b_[0];
    for(auto j=1;j<n;j++)
    {
        b[j]=b_[j]-c[j-1]*a[j]/b[j-1];
        d[j]=d[j]-d[j-1]*a[j]/b[j-1];
    }
    // calculate solution
    temp[n-1]=d[n-1]/b[n-1];
    for(int j=n-2; j>=0; j--)
        temp[j]=(d[j]-c[j]*temp[j+1])/b[j];
    return temp;
}





//if (method == PSOR){
//    // Solve matrix equations with SOR
//    int y, sor, iterMax = 10000;
//    for (sor = 0; sor < iterMax; sor++)
//    {
//        // SOR equations in here
//        // j = 0
//        {
//          y = (d[0] - c[0] * vNew[1]) / b[0];
//          vNew[0] = max(vNew[0] + omega * (y-vNew[0]), m_P);
//        }
//        // 0 < j < jMax
//        for(int j=1; j<m_J; j++)
//        {
//          y = (d[j] - a[j] * vNew[j-1] - c[j]*vNew[j+1]) / b[j];
//          vNew[j] = max(vNew[j] + omega * (y-vNew[j]), m_P);
//        }
//        // j = jMax
//        {
//          y = (d[m_J] - a[m_J] * vNew[m_J-1]) / b[m_J];
//          vNew[m_J] = max(vNew[m_J] + omega * (y-vNew[m_J]), m_P);
//        }
//        // Calculate residual
//        double error=0.;
//        error += fabs(d[0] - b[0] * vNew[0] - c[0] * vNew[1]);
//        for(int j=1; j<m_J;j++)
//          error += fabs(d[j] - a[j]*vNew[j-1] - b[j]*vNew[j] - c[j]*vNew[j+1]);
//        error += fabs(d[m_J] - a[m_J]*vNew[m_J-1] - b[m_J]*vNew[m_J]);
//        // Check for convergence and exit loop if converged
//        if(error < tol)
//          break;
//    }
//}
