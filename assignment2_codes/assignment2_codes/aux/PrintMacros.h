//
//  PrintMacros.h
//  ComputationalFinance
//
//  Created  on 24/4/21.
//


#define COMMA               << ", " <<
#define START_LINE          cout <<
#define END_LINE            << endl;

#define PRINT_DATA_LINE(a)              START_LINE (a) END_LINE
#define PRINT_2DATA_LINE(a,b)           START_LINE (a) COMMA (b) END_LINE
#define PRINT_3DATA_LINE(a,b,c)         START_LINE (a) COMMA (b) COMMA (c) END_LINE
#define PRINT_4DATA_LINE(a,b,c,d)       START_LINE (a) COMMA (b) COMMA (c) COMMA (d) END_LINE
#define PRINT_5DATA_LINE(a,b,c,d,e)     START_LINE (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) END_LINE
#define PRINT_6DATA_LINE(a,b,c,d,e,f)   START_LINE (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) END_LINE
#define PRINT_7DATA_LINE(a,b,c,d,e,f,g) START_LINE (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) END_LINE

