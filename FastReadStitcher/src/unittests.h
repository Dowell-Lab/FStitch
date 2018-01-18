#ifndef UNITTESTS_H
#define UNITTESTS_H
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cstring>
#include <cmath>

#include "main_segment.h"
#include "main_train.h"
#include "ParamWrapper.h"
#include "split.h"
#include "read.h"

typedef struct
{
    string a;
    char aCh[80];
    int aFd;
    FILE *aFile;
    
    string b;
    char bCh[80];
    int bFd;
    FILE *bFile;
} splitoutput_t;

/* genOnOff -- Splits an input training file into on and off regions to further test parameters. 
 * In the return struct:
 *  a = the "on" training example file.
 *  b = the "off" training example file.
 */
splitoutput_t *genOnOff(string infile);

/* genPosNeg -- Splits an input bedgraph into sorted positive and negative files to further test parameters.
 * In the return struct:
 *  a = the positive strand data file.
 *  b = the negative strand data file.
 */
splitoutput_t *genPosNeg(string infile);

/* checkWeightsConsistency -- Checks a given set of trained weights against a "golden standard" to determine if they are similar enough.*/
bool checkWeightsConsistency(string goldStandard, string checkFile, double errorMargin);

/* checkBedsConsistency -- Checks a given output bed file against a "golden standard" bedgraph to determine if they are similar enough.*/
bool checkBedsConsistency(string goldStandard, string checkFile);

/* cleanupSplitOutput -- Deletes any temporary files associated with a given splitoutput_t object. */
void cleanupSplitOutput(splitoutput_t *o);

void printUnitUsage(char *unitName);

#endif 
