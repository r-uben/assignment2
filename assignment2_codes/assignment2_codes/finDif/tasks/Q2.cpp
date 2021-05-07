//
//  Q2.cpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 7/5/21.
//

#include "Q2.h"
#include <vector>
#include <cmath>
#include <algorithm>

#include "CrankNicolson.h"
#include "GeneralFunctions.h"

#include <iostream>
#include <fstream>
using namespace std;

Q2::CQ2(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double Smax, int J, int I){
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
Q2::increasingS(int I, int J, int iterMax, double gap, double r)
{
    if (r == 0)
    {
        ofstream output;
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task2/amConvBondValues_I" + to_string(I) + "_J" + to_string(J) + "_iterMax" + to_string(iterMax) + "_rho1e8.csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S < iterMax; S*=gap)
        {
            CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, J, I);
            crank.amConvertibleBond(&output, PENALTY);
        }
        output.close();
    }
    else
    {
        string strR = AUX::doubleToString(r);
        ofstream output;
        output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task2/amConvBondValues_I" + to_string(I) + "_J" + to_string(J) + "_iterMax" + to_string(iterMax) + "_rho1e8r" + strR + ".csv");
        output << "F,I,J,S,V,VtoINF" << endl;
        for (double S = 2; S < iterMax; S*=gap)
        {
            CN crank(m_T, m_F, m_R, r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, J, I);
            crank.amConvertibleBond(&output, PENALTY);
        }
        output.close();
    }
}

void
Q2::variousInterestRates(vector<double> &rs, int I, int J, int iterMax, double gap)
{
    for (auto r : rs)increasingS(I, J, iterMax, gap, r);
}

void
Q2::fixedS0(double S0, double increment, int nMin, int nMax, int deg, int timesX)
{
    ofstream output;
    string strDeg = AUX::doubleToString(deg, true);
    string strSmax= AUX::doubleToString(timesX, true);
    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task2/amConvBondValues_increasing_iMax_and_jMax_deg" + strDeg + "_Smax" + strSmax + "X.csv");
    output << "F,I,J,S,V,VtoInf" << endl;
    for (int n=nMin; n<=nMax; n*=increment)
    {
        CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S0, timesX * m_X, n, n);
        crank.amConvertibleBond(&output, PENALTY, deg);
    }
}
