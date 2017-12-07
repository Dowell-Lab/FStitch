#include "main_train.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "grabTrainingExamples.h"
#include "interval_tree.h"
#include "NewtonsMethod.h"
#include "BaumWelch.h"
#include "write.h"
#include "validate.h"
using namespace std;

int run_main_train_pwrapper(ParamWrapper *p)
{
    //=================================================================
    // General Parameters

    string BedGraphFile = p->readFileName;
    string TrainingFile = p->specialFileName;
    
    string outFile = p->outFileName;
    int num_proc = p->numProcs;
    bool verbose = p->verbose;
    bool ChIP = p->chip;
    //=================================================================
    // parameters specific to Baum Welch and Newtons Method
    double max_convergence = p->maxConvergenceIters; //This was changed in the new version. ParamWrapper needs to be altered
    // to compensate!
    double convergence_threshold = p->convergenceThreshold;
    double learning_rate = p->learningRate;
    int maxSeed = p->maxSeed;

    //=================================================================
    // Misc Parameters
    string commandLine=p->commandLine;

    //=================================================================

    if (verbose){
        cout<<"reading training file                     : ";
    }
    //=================================================================
    //READ TRAINING FILE
    readTrainingFileReturn TrainReturn;
    
    if(p->strand==STRAND_BOTH)
    {
    }
    
    else if(p->strand==STRAND_POSITIVE)
    {
    }
    
    else if(p->strand==STRAND_NEGATIVE)
    {
    }

    if(p->specialFileSplit)
    {
        TrainReturn=readSplitTrainingFile(TrainingFile, p->secondSpecialFileName);
    }

    else
    {
        TrainReturn=readTrainingFile(TrainingFile);
    }

    if (TrainReturn.EXIT){
        cout<<"exiting..."<<endl;
        return 0;
    }
    map<string, interval *> TrainingIntervals = TrainReturn.result;

    if (verbose){
        cout<<"done"<<endl;
        cout<<"making interval tree                      : ";
    }
    //=================================================================
    //INTERVAL TREE FROM TRAINING FILE
    map<string,T> R = makeIntervalTree(TrainingIntervals);
    if (verbose){
        cout<<"done"<<endl;
        cout<<"reading bed graph file                    : ";
        cout<<flush;
    }
    //=================================================================
    //READ BEDGRAPH FILE
    map<string, contig*> ContigData;
    
    //This will read only positive data from the given bedgraph:
    //TODO: Implement proper support for training on both strands.
    if(p->strand==STRAND_POSITIVE || p->strand==STRAND_BOTH)
    {
        ContigData=readBedGraphFileStrand(BedGraphFile, TrainingIntervals, 1, 0);
    }
    
    else if(p->strand==STRAND_NEGATIVE)
    {
        ContigData=readBedGraphFileStrand(BedGraphFile, TrainingIntervals, 1, 1);
    }
    
    //map<string,contig *> ContigData = readBedGraphFileStrand(BedGraphFile,TrainingIntervals,1,);
    if (ContigData.empty()){
        cout<<"Set of contigs was empty. Exiting..."<<endl;
        return 0;
    }
    if (verbose){
        cout<<"done\n";
        cout<<"grabbing training data from bedgraph file : ";
        cout<<flush;

    }
    //=================================================================
    //GET DATA FROM TRAINING INTERVALS
    run_out RO = run_grabTrainingExamples(R, ContigData, ChIP);
    if (RO.EXIT){
        cout<<"exiting..."<<endl;
        return 0;
    }
    if (verbose){
        cout<<"done\n";
        cout<<"Begin parameter estimation, Newtons Method: ";
        cout<<flush;

    }
    //=================================================================
    //NEWTONS METHOD
    vector<double> W 							= learn(RO.X, RO.Y, 0, learning_rate);
    if (verbose){
        cout<<"done\n";
        cout<<"Parameter estimation, Baum-Welch          : ";
        cout<<flush;
    }
    //=================================================================
    //BAUM WELCH ALG
    BW_OUT BWO 									= runBW(ContigData, W,max_convergence, convergence_threshold,learning_rate, verbose, num_proc, maxSeed, ChIP);
    if (verbose){
        cout<<"done\n";
        cout<<"Writing learned parameters                : ";
        cout<<flush;

    }
    writeTrainingFile(outFile, BWO, learning_rate, max_convergence, convergence_threshold, ChIP, commandLine);
    if (verbose){
        cout<<"done\n";
    }
}
