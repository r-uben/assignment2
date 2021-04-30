//
//  ConvBondEDP.cpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 30/4/21.
//

#include "ConvBondEDP.h"
#include <cmath>

CB_EDP::CConvBondEDP(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma)
:CEDP(T,  r, sigma, X)
{
    m_F         = F;
    m_R         = R;
    m_kappa     = kappa;
    m_mu        = mu;
    m_C         = C;
    m_alpha     = alpha;
    m_beta      = beta;
    m_kappar    = kappa + r;
    m_alphar    = alpha + r;
}

double
CB_EDP::aFunc(double t, int j)
{
    double first_term   = -0.25 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2*(m_beta-1));
    double second_term  =  0.25 * m_kappa * ( theta(t) / m_dS - j);
    return first_term + second_term;
}

double
CB_EDP::bFunc(double t, int j)
{
    double long_term = 0.5 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2.*(m_beta-1.));
    return 1. / m_dt + 0.5 * m_r  + long_term;
}

double
CB_EDP::cFunc(double t, int j)
{
    double first_term   = -0.25 * pow(m_sigma, 2.) * pow(j, 2.*m_beta) * pow(m_dS, 2.*(m_beta-1.));
    double second_term  = -0.25 * m_kappa * ( theta(t) / m_dS - j);
    return first_term + second_term;
}

double
CB_EDP::dFunc(double t, int j, vector<double> &v)
{
    double a = aFunc(t, j);
    double b = bFunc(t, j);
    double c = cFunc(t, j);
    double d =  - ( a * v[j-1] + (b - 2 / m_dt) * v[j] + c * v[j+1] ) + m_C * exp(-m_alpha * t);
    return d;
}

// tÚ NUNCA VAS A PODER GENERAR UN OBJETO DE TIPO EDP, PERO UN PUNTERO DE TIPO EDP PUEDE APUNTAR A CUALQUIER HIJO DE EDP
// SIEMPRE QUE TÚ QUIERAS UNA FUNCIÓN QUE TRABAJE CON CLASES EDP EN GENERAL VAS A TENER QUE HACERLA DE TAL MANERA QUE SU INPUT SEA UN PUNTERO A EDP
// EN EL CONSTRUCTOR DE LA CLASE CRANK NICOLSON HAY QUE PENSAR QUÉ COSAS SON ESPECÍFICAS DE ELLA PARA QUE NOS LA PUEDA DAR LA OTRA O NO

