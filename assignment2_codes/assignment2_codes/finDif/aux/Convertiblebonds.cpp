//
//  ConvertibleBonds.cpp
//  assignment2_codes
//
//  Created on 28/4/21.
//

#include "ConvertibleBonds.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

CONV_BONDS::CConvertibleBonds(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double Smax, long I, long J)
{
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
    m_Smax      = Smax;
    m_J         = J;
    m_dS        = m_Smax / m_J;
    m_I         = I;
    m_dt        = m_T / m_I;
}

double
CONV_BONDS::A(double t)
{
    return m_R * exp(-m_kappar * (m_T - t));
}

double
CONV_BONDS::B(double t)
{
    double X_term = 2. * m_X * m_R * exp( - (m_r + 0.5 * m_kappa) * (m_T - t) ) * sinh( 0.5 * m_kappa * (m_T -t) );
    double C_term = - m_C / m_alphar * exp(m_r * t) * ( exp(-m_alphar * m_T) - exp(-m_alphar * t));
    return X_term + C_term;
}

// GENERAL VALUE FUNCTION
double
CONV_BONDS::V_Smax(double S, double t)
{
    return S * A(t) + B(t);
}

// VALUE FUNCTION WHEN S = 0
double
CONV_BONDS::V_S0(double t)
{
    double F_term = m_F * exp( -m_r * (m_T - t) );
    double C_term = m_C / ( m_r + m_alpha ) * exp( -m_alpha * m_T ) * ( exp(m_alpha * (m_T - t)) - exp(-m_r * (m_T - t)) );
    return F_term + C_term;
}

// PDE COEFFICIENTS
double
CONV_BONDS::aFunc(long i, long j)
{
    double first_term, second_term;
    double approx_t     = (i + 0.5) * m_dt;
    first_term   = -0.25 * pow(m_sigma, 2.) * pow(j, 2.* m_beta) * pow(m_dS, 2*(m_beta-1));
    second_term  =  0.25 * m_kappa * (theta(approx_t) / m_dS - j);
    return first_term + second_term;
}

double
CONV_BONDS::bFunc(long i, long j)
{
    double long_term = 0.5 * pow(m_sigma, 2.) * pow(j, 2. * m_beta) * pow(m_dS, 2. * (m_beta - 1.));
    return 1. / m_dt + 0.5 * m_r  + long_term;
}

double
CONV_BONDS::cFunc(long i, long j)
{
    double first_term, second_term;
    double approx_t     = (i + 0.5) * m_dt;
    first_term   = -0.25 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2. * (m_beta-1.));
    second_term  = -0.25 * m_kappa * (theta(approx_t) / m_dS - j);
    return first_term + second_term;
}

double
CONV_BONDS::dFunc(long i, long j, vector<double> &v)
{
    double approx_t = (i + 0.5) * m_dt;
    double a = aFunc(i, j);
    double b = bFunc(i, j);
    double c = cFunc(i, j);
    double d = - (a * v[j-1] + (b - 2. / m_dt) * v[j] + c * v[j+1]) + m_C * exp(-m_alpha * approx_t);
    return d;
}

double
CONV_BONDS::theta(double t)
{
    return (1. + m_mu) * m_X * exp(m_mu * t);
}

// LU COEFFICIENTS
double
CONV_BONDS::betaFunc(long i, long j, double prevBeta)
{
    return bFunc(i, j) - (aFunc(i, j) * cFunc(i, j-1)) / prevBeta;
}
double
CONV_BONDS::DFunc(long i, long j, double prevBeta, double d, double prevD)
{
    return d - (aFunc(i, j) * prevD) / prevBeta;
}
double
CONV_BONDS::prevV(long i, long j, vector<double> &beta, vector<double> &D, vector<double> &V)
{
    return 1. / beta[j] * (D[j] - cFunc(i, j) * V[j+1]);
}
