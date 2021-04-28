//
//  PayoffFunctions.hpp
//  ComputationalFinance
//
//  Created by Rubén Fernández Fuertes on 14/4/21.
//

#define  PAYOFF CPayoffFunctions

class CPayoffFunctions{
public:
    static double callOption(double stock_price, double strike_price);
    static double putOption(double stock_price, double strike_price);
    static double convBond(double stock_price, double R, double F);
};
