//
//  Q1.cpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 6/5/21.
//

#include "Q1.h"
#include <vector>
#include <cmath>
#include <algorithm>

#include "CrankNicolson.h"
#include "GeneralFunctions.h"

#include <iostream>
#include <fstream>
using namespace std;

Q1::CQ1(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double Smax, int J, int I){
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
    m_Smax  = Smax;
    m_J     = J;
    m_I     = I;
    m_dS    = m_Smax / J;
    m_dt    = m_T / I;
}

void
Q1::variousBetas(vector<double> &betas, int method)
{
    ofstream output;
    string strSigma = doubleToString(m_sigma);
    for (auto beta : betas)
    {
        string strBeta = doubleToString(beta);
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_beta" + strBeta + "_sigma" + strSigma + ".csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S < m_Smax; S*=1.1)
        {
            CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, beta, m_sigma, S, m_Smax, m_J, m_I);
            crank.eurConvertibleBond(&output, method);
        }
        output.close();
    }
}
void
Q1::variousKappas(vector<double> &kappas, int method)
{
    ofstream output;
    string strBeta  = doubleToString(m_beta);
    string strSigma = doubleToString(m_sigma);
    for (auto kappa : kappas)
    {
        string strKappa = doubleToString(kappa);
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_beta" + strBeta + "_kappa" + strKappa + "_sigma" + strSigma + ".csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S < m_Smax; S*=1.1)
        {
            CN crank(m_T, m_F, m_R, m_r, kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, m_J, m_I);
            crank.eurConvertibleBond(&output, method);
        }
        output.close();
    }
}

void
Q1::variousSigmas(vector<double> &sigmas, int method)
{
    ofstream output;
    string strBeta = doubleToString(m_beta);
    for (auto sigma : sigmas)
    {
        string strSigma = doubleToString(sigma);
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_beta" + strBeta + "_sigma" + strSigma + ".csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S < m_Smax; S*=1.1)
        {
            CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, sigma, S, m_Smax, m_J, m_I);
            crank.eurConvertibleBond(&output, method);
        }
        output.close();
    }
}

void
Q1::increasingS(int I, int J, int iterMax, double gap,  double r)
{
    if (r == 0)
    {
        ofstream output;
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_I" + to_string(I) + "_J" + to_string(J) + "_iterMax" + to_string(iterMax) + "_rho1e8.csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S < iterMax; S*=gap)
        {
            CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, J, I);
            crank.eurConvertibleBond(&output, THOMAS, 2);
        }
        output.close();
    }
    else
    {
        string strR = AUX::doubleToString(r);
        ofstream output;
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_I" + to_string(I) + "_J" + to_string(J) + "_iterMax" + to_string(iterMax) + "_rho1e8r" + strR + ".csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S <= iterMax; S*=gap)
        {
            CN crank(m_T, m_F, m_R, r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, J, I);
            crank.eurConvertibleBond(&output, THOMAS, 2);
        }
        output.close();
    }
}

void
Q1::fixedS0(double S0, double increment, int nMin, int nMax, int deg, int timesX)
{
    ofstream output;
    string strDeg = doubleToString(deg, true);
    string strSmax= doubleToString(timesX, true);
    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task1/eurConvBondValues_increasing_iMax_and_jMax_deg" + strDeg + "_Smax" + strSmax + "X.csv");
    output << "F,I,J,S,V,VtoInf" << endl;
    for (int n=nMin; n<=nMax; n*=increment)
    {
        CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S0, timesX * m_X, n, n);
        crank.eurConvertibleBond(&output, THOMAS, deg);
    }
}

string
Q1::doubleToString(double value, bool integer)
{
    string str = to_string(value);
    if (value < 10)
        str = str.erase(1,1);
    if (value > 10 && value < 100)
        str = str.erase(2,2);
    if (integer == true)
        str = str.erase(1);
    else
        str = str.erase(5);
    return str;
}
