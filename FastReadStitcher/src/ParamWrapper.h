#ifndef PARAMWRAPPER_H
#define PARAMWRAPPER_H

#include <iostream>
#include <cstdio>
#include <fstream>
#include <unistd.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <cstdlib>
using namespace std;

/* ParamWrapper.h -- Header for an updated parameter parsing and storage class.
 *
 * It implements roughly the same same interface as paramWrapper, but in a way that is
 * significantly more user friendly and less awful.
 */

#define STRAND_POSITIVE 1
#define STRAND_NEGATIVE 2
#define STRAND_BOTH 3
//Certain versions of FStitch support automatic determination of which
//strand a given set of parameters was trained with:
#define STRAND_UNSPECIFIED 4

#define REPORT_ON 1
#define REPORT_OFF 2
#define REPORT_BOTH 3

class ParamWrapper
{
public:
    string command;
    
    bool train, segment, eRNA, exit;
    map<string, string> params;
    
    string outFileName;
    //This represents the annotation file for train and
    //the the weights file for segment.
    string specialFileName;
    //This represents the raw read data:
    string readFileName;
    
    bool specialFileSplit;
    string secondSpecialFileName;
    
    //This string represents the command line passed to the program. It will be used
    //when writing output files so that results are easier to reproduce.
    string commandLine;
    
    int strand;
    int report;
    
    int numProcs;
    
    /* Provide some additional values related to the learning process: */
    int maxConvergenceIters;
    double convergenceThreshold;
    double learningRate;
    double regularization;
    int maxSeed;
    
    bool chip;
    bool verbose;
    
    ParamWrapper(int argc, char **argv);
    
    /* This generic constructor is used by unit tests to manipulate various options on the fly:*/
    ParamWrapper();
    void printUsage();
    void dumpValues();
};
#endif
