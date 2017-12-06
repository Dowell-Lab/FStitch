//============================================================================
// Name        : main.cpp
// Author      : Joey Azofeifa
// Version     : 1.1
// Description : Main file for running FastReadStitcher Package
//============================================================================
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include "read_in_parameters.h"
#include "main_train.h"
#include "main_segment.h"
#include "ParamWrapper.h"
using namespace std;

int main(int argc, char* argv[]) {
    ParamWrapper *p=new ParamWrapper(argc, argv);
    if(p->exit)
    {
        cout<<"exiting..."<<endl;
        delete p;
        return 0;
    }

    if(p->verbose)
    {
        cout<<"Verbose specified. Printing usage for some reason:"<<endl;
        p->printUsage();
    }

    if(p->train)
    {
        run_main_train_pwrapper(p);
    }

    else if(p->segment)
    {
        run_main_segment_pwrapper(p);
    }

    else if(p->eRNA)
    {
        cout<<"need to finish segment code"<<endl;
    }

    delete p;
    return EXIT_SUCCESS;
}
