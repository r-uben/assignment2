//
//  EDP.h
//  assignment2_codes
//
//  Created on 30/4/21.
//

#include <vector>
#include <fstream>
using namespace std;
d
class CEDP
{
public:
    CEDP(double T, double r, double sigma, double X);
    /// To do:  definir estructura para la opción
    virtual double aFunc(double t, int j) = 0;
    virtual double bFunc(double t, int j) = 0;
    virtual double cFunc(double t, int j) = 0;
    virtual double dFunc(double t, double j, vector<double>& v) = 0;
    
private:
    double m_T;
    double m_r;
    double m_sigma;
    double m_X;
};
