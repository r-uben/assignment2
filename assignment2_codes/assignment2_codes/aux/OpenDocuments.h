//
//  OpenDocuments.hpp
//  assignment2_codes
//
//  Created by Rubén Fernández Fuertes on 30/4/21.
//

#include <iostream>
#include <fstream>
using namespace std;

void OpenCSVFile(ofstream *output, string name){
    // open up a file stream to write data
    (*output).open("/Users/rubenexojo/Library/Mobile Documents/com~apple~CloudDocs/MSc Mathematical Finance - Manchester/subjects/II_semester/MATH60082_computational_finance/c++/assignment2/data/" + name + ".csv");
    // check if the file is opened
    if (!(*output).is_open()){
        cout << " File not opened" << endl;
        // stop here
        throw;
    }
}
