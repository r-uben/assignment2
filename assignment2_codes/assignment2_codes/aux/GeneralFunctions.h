//
//  GeneralFunctions.h
//  ComputationalFinance
//
//  Created on 14/4/21.
//

#include <vector>

#define  AUX CGeneralFunctions

using namespace std;

class CGeneralFunctions{
public:
    static double discountFactor(double interest_rate, double time);
    static double maxFunc(double x, double y);
    static double minFunc(double x, double y);
    static double lagInterp(double S, vector<double> p1, vector<double> p2);
};
