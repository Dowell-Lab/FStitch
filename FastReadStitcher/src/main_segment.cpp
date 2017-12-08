#include "main_segment.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "interval_tree.h"
#include "write.h"
#include "viterbi.h"
#include "validate.h"

#include <string>
using namespace std;
bool isNeg(string OUT){
	for (int i = 0; i < OUT.size(); i++){
		for (int j = i; j < OUT.size(); j++){
			if (OUT.substr(i,OUT.size()-j)=="neg" or OUT.substr(i,OUT.size()-j) =="-"){
				return 1;
			}
		}
	}
	return 0;
}
bool isPos(string OUT){
	for (int i = 0; i < OUT.size(); i++){
		for (int j = i; j < OUT.size(); j++){
			if (OUT.substr(i,OUT.size()-j)=="pos" or OUT.substr(i,OUT.size()-j) =="+"){
				return 1;
			}
		}
	}
	return 0;
}

int run_main_segment_pwrapper(ParamWrapper *p)
{
    //=================================================================
    // General Parameters
    string BedGraphFile 		= p->readFileName;
    string TrainingOutFile 		= p->specialFileName;
    string outFile 				= p->outFileName;
    //Todo: add functionality for this!
    //string refFile=p->params["-r"];
    string refFile;
    //string refFile 				= PT.params["-r"];
    //This will be handled differently:
    string strand;

    //==================================================================
    // Other Parameters
    string commandLine=p->commandLine;

    if(p->strand==STRAND_UNSPECIFIED)
    {
        if (isPos(BedGraphFile)){
            strand = "+";
        }else if(isNeg(BedGraphFile)){
            strand = "-";
        }else{
            strand = ".";
        }
    }

    else if(p->strand==STRAND_BOTH)
    {
        strand=".";
    }

    else if(p->strand==STRAND_POSITIVE)
    {
        strand="+";
    }

    else if(p->strand==STRAND_NEGATIVE)
    {
        strand="-";
    }

    int num_proc = p->numProcs;

    if (strand.empty()){
        if (isPos(BedGraphFile)){
            strand = "+";
        }else if(isNeg(BedGraphFile)){
            strand = "-";
        }else{
            strand = ".";
        }
    }
    bool verbose=p->verbose;
    //=================================================================
    //Read in FStich Training out file
    if (verbose){
        cout<<"Reading in FStitch Training Out File: ";
    }
    RTOF RTOF_params 					= readTrainingOutFile(TrainingOutFile);

    if (RTOF_params.EXIT){
        cout<<"exiting..."<<endl;
        return 0;
    }
    if (verbose){
        cout<<"done"<<endl;
        cout<<flush;
    }
    if (verbose){
        cout<<"Training Set from ChIP Data         : "<<(bool(RTOF_params.ChIP)==1)<<endl;
    }
    //=================================================================
    //Read in BedGraph File by Chromosome
    if (verbose){
        cout<<"Reading in BedGraph File            : ";
        cout<<flush;
    }
    map<string,contig *> ContigData 	= readBedGraphFileAll(BedGraphFile,num_proc);
    if (ContigData.empty()){
        cout<<"exiting..."<<endl;
        return 0;
    }

    if (verbose){
        cout<<"done"<<endl;
        cout<<flush;
    }
    //=================================================================
    //Run Viterbi
    if (verbose){
        cout<<"Running Viterbi                     : ";
        cout<<flush;
    }
    map<string, state*> results 	= runViterbi(ContigData, RTOF_params.W, RTOF_params.A,num_proc, RTOF_params.ChIP);
    if (verbose){
        cout<<"done"<<endl;
        cout<<flush;
    }
    //=================================================================
    //Write Viterbi Paths
    if (verbose){
        cout<<"Writing to IGV                      : ";
        cout<<flush;
    }
    writeViterbiPaths(outFile, results, refFile, strand, RTOF_params.commandLine, commandLine);
    if (verbose){
        cout<<"done"<<endl;
    }

    return 1;
}
