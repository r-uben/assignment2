//
//  GeneralFunctions.cpp
//  ComputationalFinance
//
//  Created on 14/4/21.
//

#include "GeneralFunctions.h"
#define  AUX CGeneralFunctions

#include <cmath>
#include <vector>
using namespace std;

double
AUX::discountFactor(double interest_rate, double time)
{
    return exp(- interest_rate * time);
}
double
AUX::maxFunc(double x, double y)
{
    if (x <= y)
        return y;
    else
        return x;
}
double
AUX::minFunc(double x, double y)
{
    if (x <= y)
        return x;
    else
        return y;
}
double
AUX::lagInterp(double S, vector<double> p1, vector<double> p2)
{
    double approxJ        = (S - p2[0]) / (p1[0] - p2[0]) * p1[1];
    double approxJplus1   = (S - p1[0]) / (p2[0] - p1[0]) * p2[1];
    return approxJ + approxJplus1;
}
