//
//  GridConstructor.h
//  ComputationalFinance
//
//  Created by Rubén Fernández Fuertes on 15/4/21.
//

#include <vector>
#define  GRID CGridConstructor

using namespace std;

class CGridConstructor{
public:
    static vector<double>   setupStockPrices(double dS, int J);
    static double           boundaryConditions(int j);
    static double           lagrangeInterpolation(const vector<double>& y,const vector<double>& x, double x0, unsigned int n);
};
