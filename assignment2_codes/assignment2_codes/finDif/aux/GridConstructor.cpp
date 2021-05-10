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
GRID::setupStockPrices(double dS, long J){
    vector<double> S(J+1);
    for (int j=0; j<=J; j++){
        S[j]  = j * dS;
    }
    return S;
}

// A generic llgrange interpolation function
double
GRID::lagrangeInterpolation(const vector<double>& y,const vector<double>& x,double x0, unsigned int n)
    {
        if(x.size()<n)return lagrangeInterpolation(y,x,x0,x.size());
        if(n==0)throw;
        int nHalf = n/2;
        int jStar;
        double dx=x[1]-x[0];
        if(n%2==0)
            jStar = int((x0 - x[0])/dx) -(nHalf-1);
        else
            jStar = int((x0 - x[0])/dx+0.5)-(nHalf);
        jStar= max(0,jStar);
        jStar= min(int(x.size()-n),jStar);
        if(n==1)
            return y[jStar];
        double temp = 0.;
        for(unsigned int i=jStar;i<jStar+n;i++){
            double  int_temp;
            int_temp = y[i];
            for(unsigned int j=jStar;j<jStar+n;j++){
                if(j==i){continue;}
                int_temp *= ( x0 - x[j] )/( x[i] - x[j] );
            }
            temp += int_temp;
        }
        // end of interpolate
        return temp;
    }
