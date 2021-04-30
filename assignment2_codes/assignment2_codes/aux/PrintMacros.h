//
//  PrintMacros.h
//  ComputationalFinance
//
//  Created  on 24/4/21.
//

#include <iostream>
#include <fstream>
using namespace std;

#define COMMA               << ", "
#define START_LINE          cout <<
#define END_LINE            << endl;

/// CONSOLE
// RECURSION FOR PRINTING RESULTS IN THE CONSOLE: WE FIRSTLY DECLARE THE "LAST STATE" IN THE
// RECURSION, I.E., WHEN THERE'S ONLY ONE ARGUMENT LEFT TO PRINTING.
// THE SECOND TEMPLATE DOES THE REST.

template<typename T>                    // Type is resolved in compile time
void PRINT_DATA_LINE(T t)
{
    START_LINE t END_LINE;
}

template<typename T, typename... ARGS>
void PRINT_DATA_LINE(T t, ARGS... args) // Take the first arguments
{
    START_LINE t COMMA;
    PRINT_DATA_LINE(args...);                 // Recursion with a fewer argument
}

/// CSV
// (SAME PROCEDURE THAN BEFORE)
template<typename T>
void DATA_LINE(ofstream *output, T t)
{
    *output << t << endl;
}

template<typename T, typename... ARGS>
void DATA_LINE(ofstream *output, T t, ARGS... args)
{
    *output << t << ", ";
    DATA_LINE(output, args...);
}

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


