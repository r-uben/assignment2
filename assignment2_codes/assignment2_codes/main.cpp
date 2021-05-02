
#include <cmath>
#include <vector>
#include <algorithm>

#include "CrankNicolson.h"

#include <iostream>
#include <fstream>
using namespace std;

using namespace std;
/* Template code for the Crank Nicolson Finite Difference
 */
int main(){
    // PARAMETERS

    double T = 3., F = 25., R = 2, r = 0.0117, kappa = 0.0833333333333, mu = 0.0053, X = 17.38, C = 0.205, alpha = 0.02, beta = 0.808, sigma = 0.66;
    // double T = 2., F = 310., R = 4, r = 0.021, kappa = 0.125, mu = 0.0198, X = 77.36, C = 3.26, alpha = 0.02, beta = 0.528, sigma = 2.91;
    double Smax = R*F;
    double I = 300;
    double J = 300;

    ofstream output;
    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/eurConvBondValues.csv");

    output << "F,S,V,VtoINF" << endl;

    cout << "################"<< endl;
//    CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 75, Smax, jMax, iMax);
//    crank.convertibleBond(&output, LU);
    for (double S = 2; S < Smax; S*=1.1)
    {
        CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta=0.808, sigma=0.66, S, Smax, J, I);
        crank.convertibleBond(&output, THOMAS);
    }
    output.close();

    return 0;
}
