//
//  PayoffFunctions.cpp
//  ComputationalFinance
//
//  Created on 14/4/21.
//

#include "PayoffFunctions.h"
#include "GeneralFunctions.h"

// CALL OPTION PAYOFF
double
PAYOFF::callOption(double stock_price, double strike_price){
    return AUX::maxFunc(stock_price - strike_price, 0.);
}

// PUT OPTION PAYOFF
double
PAYOFF::putOption(double stock_price, double strike_price){
    return AUX::maxFunc(strike_price - stock_price, 0.);
}

// CONVERTIBLE BOND
double
PAYOFF::convBond(double stock_price, double R, double F)
{
    return AUX::maxFunc(F, R * stock_price);
}
