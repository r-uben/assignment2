//
//  ConvertibleBonds.hpp
//  assignment2_codes
//
//  Created on 28/4/21.
//

#define CONV_BONDS CConvertibleBonds

class CConvertibleBonds
{
public:
    CConvertibleBonds(double T, double F, double R, double r, double kappa, double mu, double X, double C, double alpha, double beta, double sigma);
    
    double V(double S, double t);
    double V_S0(double t);
    // AUXILIAR FUNCTIONS
    double A(double t);
    double B(double t);
    
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
    // AUXILIAR PARAMETERS
    double m_kappar;
    double m_alphar;
};
