//
//  PrintMacros.h
//  ComputationalFinance
//
//  Created  on 24/4/21.
//

#include <iostream>
#include <fstream>
using namespace std;

#define COMMA               << ", " <<
#define START_LINE          cout <<
#define END_LINE            << endl;

// CONSOLE
#define PRINT_DATA_LINE(a)              START_LINE (a) END_LINE
#define PRINT_2DATA_LINE(a,b)           START_LINE (a) COMMA (b) END_LINE
#define PRINT_3DATA_LINE(a,b,c)         START_LINE (a) COMMA (b) COMMA (c) END_LINE
#define PRINT_4DATA_LINE(a,b,c,d)       START_LINE (a) COMMA (b) COMMA (c) COMMA (d) END_LINE
#define PRINT_5DATA_LINE(a,b,c,d,e)     START_LINE (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) END_LINE
#define PRINT_6DATA_LINE(a,b,c,d,e,f)   START_LINE (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) END_LINE
#define PRINT_7DATA_LINE(a,b,c,d,e,f,g) START_LINE (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) END_LINE

// CSV
#define OUT                                                     output <<
#define OUT_DATA_LINE_2(out, a, b)                              out (a) COMMA (b) END_LINE
#define OUT_DATA_LINE_3(out, a, b, c)                           out (a) COMMA (b) COMMA (c) END_LINE
#define DATA_LINE_2(a, b)                                       OUT_DATA_LINE_2(OUT, a, b)
#define DATA_LINE_3(a, b, c)                                    OUT_DATA_LINE_3(OUT, a, b, c)
#define DATA_LINE_4(a, b, c, d)                                 OUT (a) COMMA (b) COMMA (c) COMMA (d) END_LINE
#define DATA_LINE_5(a, b, c, d, e)                              OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) END_LINE
#define DATA_LINE_6(a, b, c, d, e, f)                           OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) END_LINE
#define DATA_LINE_7(a, b, c, d, e, f, g)                        OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) END_LINE
#define DATA_LINE_8(a, b, c, d, e, f, g, h)                     OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) COMMA (h) END_LINE
#define DATA_LINE_9(a, b, c, d, e, f, g, h, i)                  OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) COMMA (h) COMMA (i) END_LINE
#define DATA_LINE_10(a, b, c, d, e, f, g, h, i, j)              OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) COMMA (h) COMMA (i) COMMA (j) END_LINE
#define DATA_LINE_11(a, b, c, d, e, f, g, h, i, j, k)           OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) COMMA (h) COMMA (i) COMMA (j) COMMA (k) END_LINE
#define DATA_LINE_12(a, b, c, d, e, f, g, h, i, j, k, l)        OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) COMMA (h) COMMA (i) COMMA (j) COMMA (k) COMMA (l) END_LINE
#define DATA_LINE_13(a, b, c, d, e, f, g, h, i, j, k, l, m)     OUT (a) COMMA (b) COMMA (c) COMMA (d) COMMA (e) COMMA (f) COMMA (g) COMMA (h) COMMA (i) COMMA (j) COMMA (k) COMMA (l) COMMA (m) END_LINE


// TEX

#define AMP                             " & "  <<
#define NEW_LINE                        " \\"  << "\\" << endl;
#define START_TABULAR                   output << "\\" << "begin{tabular}{cccc}" << endl;
#define END_TABULAR                     output << "\\" << "end{tabular}" << endl;
#define WRITE_HEADER(a, b, c, d)        OUT (a)  << AMP (b) << AMP (c) << AMP (d)  << " \\" << "\\" << "\\hline" << endl;
#define WRITE_TABLE_LINE(a, b, c, d)    OUT (a)  << AMP (b) << AMP (c) << AMP (d)  << NEW_LINE
#define START_EQUATION                  output << "\\begin{equation}" << endl;
#define WRITE_OPTION_VALUE(t, V)        output << "\\boxed{ V(r_0, t = " << (t) << ", T) =" <<  (V) << "}"  << endl;
#define END_EQUATION                    output << "\\end{equation}" << endl;

#define OVER_WRITE true
#define NOT_OVER_WRITE false

void OpenCSVFile(ofstream *output, string name, bool over_write){
    // open up a file stream to write data
    if (over_write == true)
        (*output).open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/" + name + ".csv");
    if (over_write ==
        false)
        (*output).open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/" + name + ".csv", fstream::app);
    // check if the file is opened
    if (!(*output).is_open()){
        cout << " File not opened" << endl;
        // stop here
        throw;
    }
}
