#ifndef write_H
#define write_H
#include <string>
#include <map>

/* These are to support setting extended attributes on files to represent the input commandline. */
#include <sys/types.h>
#include <attr/xattr.h>
#include <cstring>

#include "BaumWelch.h"
#include "viterbi.h"
void writeTrainingFile(string, BW_OUT,double, double, double, bool, string);
void writeViterbiPaths(string, map<string, state*>, string, string, string);
#endif