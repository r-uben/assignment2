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

// To keep results
#define COMMA << "," <<
#define END_LINE << endl;

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
            CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, I, J);
            crank.amConvertibleBond_penalty(&output,2);
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
            CN crank(m_T, m_F, m_R, r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S, m_Smax, I, J);
            crank.amConvertibleBond_penalty(&output, 2);
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
Q2::variousV_fixedS0(double S0, double t0, double increment, int nMin, int nMax, int deg, int timesX)
{
    double Smax = timesX * m_X;
    ofstream output;
    // Convert the integers into strings to name the document
    string strDeg = to_string(deg);
    string strSmax= to_string(timesX);
    // Open the document in order to keep results
    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task2/amConvBondValues_increasing_iMax_and_jMax_deg" + strDeg + "_Smax" + strSmax + "X.csv");
    // First Line of the .csv file to further get columns of data
    output << "F,I,J,S,V,VtoInf" << endl;
    // Iteration on the size of the grid. I and J are conveniently chosen in order to have S0 and t0 on the grid
    // or close to the grid.
    for (int n=nMin; n<=nMax; n*=increment)
    {
        // Produce and keep results
        CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S0, Smax, n * ceil(Smax / S0), n * ceil(m_T / t0));
        crank.amConvertibleBond_penalty(&output, deg);
    }
}

void
Q2::V_fixedS0(double S0, double t0, int deg, int timesX, int I, int J)
{
    // Convert the integers into strings to name the document
    ofstream output;
    // Open the document in order to keep results
    string strDeg = to_string(deg);
    string strSmax= to_string(timesX);
    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/task2/amConvBondValue_" + strSmax + "X_deg" + strDeg + "_I" + to_string(I) + "_J" + to_string(J) + "_timing.csv");
    // First Line of the .csv file to further get columns of data
    output << "I,J,V,time" << endl;
    // Keep the time of inisiation of the Crank Nicolson method
    auto start = START_TIME;
    // Produce results
    CN crank(m_T, m_F, m_R, m_r, m_kappa, m_mu, m_X, m_C, m_alpha, m_beta, m_sigma, S0, timesX * m_X, I, J);
    crank.amConvertibleBond_penalty(&output, deg, DONT_SAVE);
    double optionValue = crank.GetV();
    // Keep the time of the end of the Crank Nicolson method
    auto end = END_TIME;
    // Duration: end - start
    DURATION<float> duration = (end - start);
    output << I COMMA J COMMA optionValue COMMA duration.count() END_LINE
    cout << "Our result has been obtained in " << duration.count() << "s" << endl;
}

