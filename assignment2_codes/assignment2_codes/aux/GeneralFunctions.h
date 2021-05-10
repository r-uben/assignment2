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
    static string doubleToString(double value);
    static void   extrap(double S1, double S2, double p);
};

// To calculate the time
#include <chrono>
#define  CHRONO   std::chrono
#define  SET_TIME CHRONO::system_clock::now()
#define  START_TIME CHRONO::system_clock::now()
#define  END_TIME CHRONO::system_clock::now()
#define  DURATION CHRONO::duration
#define  MILLI    std::milli
