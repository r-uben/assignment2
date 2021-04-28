//
//  main.cpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 28/4/21.
//

#include <iostream>
#include "CrankNicolson.h"

int main()
{
    // PARAMETERS
    double T = 3., F = 25., R = 2, r = 0.0117, kappa = 0.0833333333333, mu = 0.0053, X = 17.38, C = 0.205, alpha = 0.02, beta = 0.808, sigma = 0.66;
    double S0 = 17.38, Smax = 2*F*R;
    double I = 1000;
    double J = 1000;
    
    CN crank(T, F, R, r, kappa, mu, X, C, alpha, 1, 0.381, S0, Smax, J, I);
    crank.convertibleBond(LU);
    return 0;
}
