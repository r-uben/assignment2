//
//  main.cpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 28/4/21.
//

#include "CrankNicolson.h" 

#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    // PARAMETERS
    // mine double T = 3., F = 25., R = 2, r = 0.0117, kappa = 0.0833333333333, mu = 0.0053, X = 17.38, C = 0.205, alpha = 0.02, beta = 0.808, sigma = 0.66;
    double T = 2., F = 310., R = 4, r = 0.021, kappa = 0.125, mu = 0.0198, X = 77.36, C = 3.26, alpha = 0.02, beta = 0.528, sigma = 2.91;
    double Smax = 3*F*R;
    double I = 100;
    
    double J = 100;
    
    ofstream output;
    output.open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/eurConvBondValues.csv");
    
    output << "F,S,V" << endl;

    CN crank(T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, 75, Smax, J, I);
    crank.convertibleBond(LU);

    output.close();
    
    return 0;
}
