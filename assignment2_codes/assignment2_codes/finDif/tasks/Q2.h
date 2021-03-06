//
//  Q2.hpp
//  assignment2_codes
//
//  Created on 7/5/21.
//

#define Q2 CQ2
#include <vector>
using namespace std;

class CQ2
{
public:
    CQ2(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma, double Smax, int J, int I);
    // Main functions
    void increasingS(int I, int J, int iterMax, double gap = 1.1, double r = 0);
    void variousInterestRates(vector<double> &rs, int I, int J, int iterMax, double gap = 1.1);
    void V_fixedS0(double S0, double t0, int deg, int timesX, int I, int J);
    void variousV_fixedS0(double S0, double t0, double increment, int nMin, int nMax, int deg, int timesX);
private:
    //
    string doubleToString(double value, bool integer = false);
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
