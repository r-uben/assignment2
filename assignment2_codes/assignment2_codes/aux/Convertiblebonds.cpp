//
//  ConvertibleBonds.cpp
//  assignment2_codes
//
//  Created on 28/4/21.
//

#include "ConvertibleBonds.h"
#include <cmath>
#include <algorithm>

CONV_BONDS::CConvertibleBonds(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma)
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
}

double
CONV_BONDS::A(double t)
{
    return m_R * exp(- m_kappar * (m_T - t));
}

double
CONV_BONDS::B(double t)
{
    double X_term = 2 * m_X * m_R * exp( -(m_r - m_kappar) * (m_T - t) ) * sinh( 0.5 * m_kappar * (m_T -t) );
    double C_term = - m_C / m_alphar * exp(m_r * t) * ( exp(-m_alphar * m_T) - exp(-m_alphar * t));
    return X_term + C_term;
}

double
CONV_BONDS::V(double S, double t)
{
    return S * A(t) + B(t);
}
