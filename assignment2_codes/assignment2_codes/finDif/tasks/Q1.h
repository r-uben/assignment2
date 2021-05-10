//
//  Q1.hpp
//  assignment2_codes
//
//  Created on 6/5/21.
//

#define Q1 CQ1
#include <vector>
using namespace std;

class CQ1
{
public:
    CQ1(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double Smax, int J, int I);
    // Main functions
    void variousBetas(vector<double> &betas, int method);
    void variousKappas(vector<double> &kappas, int method);
    void variousSigmas(vector<double> &sigmas, int method);
    void increasingS(int I, int J, int iterMax, double gap = 1.1,  double r = 0);
    void V_fixedS0(double S0, int deg, int timesX, int I, int J);
    void variousV_fixedS0(double S0, double increment, int nMin, int nMax, int deg, int timesX);
private:
    // Useful Functions
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
    double m_Smax;
    double m_dS;
    double m_dt;
    // Grid Parameters
    int    m_J;
    int    m_I;
    int    m_jStar;
};

// To calculate the time
#include <chrono>
#define  CHRONO   std::chrono
#define  SET_TIME CHRONO::system_clock::now()
#define  START_TIME CHRONO::system_clock::now()
#define  END_TIME CHRONO::system_clock::now()
#define  DURATION CHRONO::duration
#define  MILLI    std::milli
