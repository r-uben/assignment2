//
//  CrankNicolson.h
//  assignment2_codes
//
//  Created on 28/4/21.
//

#include <vector>
using namespace std;

#define CN  CCrankNicolson
#define LU  true
#define SOR false

class CCrankNicolson
{
public:
    CCrankNicolson(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double S0, double Smax, int J, int I);
    // Main functions
    void convertibleBond(bool lu, double tol = 1.e-8, double omega = 1.);
private:
    /// Other Functions
    // PDE coefficients
    double aFunc(int j);
    double bFunc(int j);
    double cFunc(int j);
    double dFunc(int j, vector<double> &v);
    // LU coefficients
    double betaFunc(int j, double prevBeta);
    double DFunc(int j, double prevBeta, double d, double prevD);
    double prevV(int j, vector<double>& beta, vector<double>& D, vector<double>& V);
    // Interior grid points
    double interiorPoints(vector<double> &v, int j);
    double approxPrice(vector<double> &v, vector<double> &s);
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
    // AUXILIAR PARAMETERS
    double m_kappar;
    double m_alphar;
    // Local Variables (dS, dt)
    double m_S0;
    double m_Smax;
    double m_dS;
    double m_dt;
    // Grid Parameters
    int    m_J;
    int    m_I;
    int    m_jStar;
};