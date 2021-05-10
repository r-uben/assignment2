//
//  GeneralFunctions.cpp
//  ComputationalFinance
//
//  Created on 14/4/21.
//

#include "GeneralFunctions.h"
#define  AUX CGeneralFunctions
#define COMMA << "," <<
#define END_LINE << endl;


#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

double
AUX::discountFactor(double interest_rate, double time)
{
    return exp(-interest_rate * time);
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
string
AUX::doubleToString(double value)
{
    string str = to_string(value);
    if (value < 10)
        // erase the comma of a number lower than 10
        str = str.erase(1,1);
    if (value > 10 && value < 100)
        // erase the comma of a number greater than 10
        str = str.erase(2,2);
    // 4-5 decimals
    str = str.erase(5);
    return str;
}

void
AUX::extrap(double S1, double S2, double p)
{
    auto start = START_TIME;
    double extrapValue = (pow(2,p)*S2 - S1) / (pow(2,p) - 1);
    auto end   = END_TIME;
    DURATION<float> duration = (end - start);
    cout << extrapValue COMMA duration.count() END_LINE;
}

