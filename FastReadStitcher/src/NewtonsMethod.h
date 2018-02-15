#ifndef NEWTONSMETHOD_H_
#define NEWTONSMETHOD_H_
#include <vector>
vector<double> learn(vector< vector<double> >, vector<int> , bool, double);
vector<double> ParameterizedNewtonsMethod(vector< vector<double> > X, vector<int> Y, bool verbose, double alpha, int numIters, double thresh);
double g(vector<double> , vector<double> );

#endif
