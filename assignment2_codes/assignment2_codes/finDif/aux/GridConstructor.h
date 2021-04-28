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
};
