//
//  ConvertibleBonds.hpp
//  assignment2_codes
//
//  Created on 28/4/21.
//

#define CONV_BONDS CConvertibleBonds

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class CConvertibleBonds
{
public:
    CConvertibleBonds(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double Smax, long I, long J);
    
    double V_Smax(double S, double t);
    double V_S0(double t);
    // AUXILIAR FUNCTIONS
    double A(double t);
    double B(double t);
    // COEFFICIENTS FOR THE PDE
    double aFunc(long i, long j);
    double bFunc(long i, long j);
    double cFunc(long i, long j);
    double dFunc(long i, long j, vector<double> &v);
    double theta(double t);
    // LU coefficients
    double betaFunc(long i, long j, double prevBeta);
    double DFunc(long i, long j, double prevBeta, double d, double prevD);
    double prevV(long i, long j, vector<double>& beta, vector<double>& D, vector<double>& V);
private:
    // PARAMETERS
    double m_T;
    double m_F;
    double m_R;
    double m_r;
    double m_kappa;
    double m_mu;
    double m_X;
    double m_C;
    double m_alpha;
    double m_beta;
    double m_sigma;
    // PDE PARAMETERS
    double m_Smax;
    double m_J;
    double m_I;
    double m_dS;
    double m_dt;
    // AUXILIAR PARAMETERS
    double m_kappar;
    double m_alphar;
    
};
