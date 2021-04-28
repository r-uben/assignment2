//
//  GridConstructor.cpp
//  ComputationalFinance
//
//  Created by Rubén Fernández Fuertes on 15/4/21.
//

#include "GridConstructor.h"
#include <vector>
using namespace std;

vector<double>
CGridConstructor::setupStockPrices(double dS, int J){
    vector<double> S(J+1);
    for (int j=0; j<=J; j++){
        S[j]  = j * dS;
    }
    return S;
}
