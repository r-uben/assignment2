//
//  ConvBondEDP.hpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 30/4/21.
//

#include "EDP.h"
#include <vector>
#include <fstream>
#include <cmath>

#define CB_EDP CConvBondEDP

using namespace std;

class CConvBondEDP:public CEDP{
public:
    CConvBondEDP(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma);
    double aFunc(double t, int j);
    double bFunc(double t, int j);
    double cFunc(double t, int j);
    double dFunc(double t, int j, vector<double>& v);
    
private:
    // PARAMETERS
    double m_F;
    double m_R;
    double m_kappa;
    double m_mu;
    double m_C;
    double m_alpha;
    double m_beta;
    double m_kappar;
    double m_alphar;
};
