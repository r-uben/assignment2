//
//  CrankNicolson.h
//  assignment2_codes
//
//  Created on 28/4/21.
//

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#define CN  CCrankNicolson
#define LU      1
#define SOR     2
#define THOMAS  3

#define PSOR                1
#define POLICY_ITERATION    2
#define PENALTY             3

#define EUROPEAN false
#define AMERICAN true

#define SAVE        true
#define DONT_SAVE   false

class CCrankNicolson
{
public:
    CCrankNicolson(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double S0, double Smax, long I, long J);
    // Main functions
    void eurConvertibleBond(ofstream *output, int method = 0, int degree = 2, double tol = 1.e-2, double omega = 1.2);
    void amConvertibleBond_penalty(ofstream *output, int degree=2, bool saveData = true, double tol = 1.e-10);
    // Get Option Value
    inline double GetV() {return m_optionValue;};
private:
    // Useful Functions
    double approxPrice(vector<double> &v, vector<double> &s);
    vector<double> thomasSolve(const vector<double> &a,const vector<double> &b_,const vector<double> &c, vector<double> &d);
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
    // AMERICAN PARAMETERS
    double m_P;
    double m_rho;
    double m_iterMax;
    double m_t0;
    // AUXILIAR PARAMETERS
    double m_kappar;
    double m_alphar;
    // Local Variables (dS, dt)
    double m_S0;
    double m_Smax;
    double m_dS;
    double m_dt;
    // Grid Parameters
    long   m_J;
    long   m_I;
    int    m_jStar;
    // Option Value
    double m_optionValue;
};
