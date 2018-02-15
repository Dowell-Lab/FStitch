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
    string BedGraphNegFile=p->secondReadFileName;
    bool BedGraphSplit=p->readFileSplit;
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
    string strand;
    int determinedStrand;

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
    
    if(p->strand==STRAND_UNSPECIFIED)
    {
        cout<<"NOTE: segmented strand not specified. Attempting to determine strand from reads histogram file..."<<endl;
        if (isPos(BedGraphFile)){
            strand = "+";
            cout<<"Input bedgraph file determined to be positive based on file name."<<endl;
        }else if(isNeg(BedGraphFile)){
            strand = "-";
            cout<<"Input bedgraph file determined to be negative based on file name"<<endl;
        }else{
            /* Attempt to determine the bed file type based on its contents: */
            determinedStrand=checkBedFileType(BedGraphFile);
            
            if(determinedStrand==STRAND_POSITIVE)
            {
                strand="+";
                cout<<"Input bedgraph file determined to be positive based on file contents."<<endl;
            }
            
            else if(determinedStrand==STRAND_NEGATIVE)
            {
                strand="-";
                cout<<"Input bedgraph file determined to be negative based on file contents."<<endl;
            }
            
            else if(determinedStrand==STRAND_BOTH)
            {
                strand=".";
                cout<<"Input bedgraph file determined to contain both positive and negative strand data"<<endl;
                cout<<"based on file contents."<<endl;
            }
            
            else
            {
                cout<<"Input bedgraph file is poorly formatted or unreadable. Exiting..."<<endl;
                return 0;
            }
            strand = ".";
            cout<<"Input bedgraph file determined to be both positive and negative."<<endl;
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
    
    if(p->specialFileSplit)
    {
        TrainReturn=readSplitTrainingFile(TrainingFile, p->secondSpecialFileName);
    }

    else
    {
        TrainReturn=readTrainingFile(TrainingFile);
    }

    if (TrainReturn.EXIT){
        cout<<"TrainReturn.EXIT set. exiting..."<<endl;
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
    /*
    if(p->strand==STRAND_POSITIVE || p->strand==STRAND_BOTH)
    {
        if(BedGraphSplit)
        {
            ContigData=readSplitBedGraphFileStrand(BedGraphFile, BedGraphNegFile, TrainingIntervals, 1, 0);
        }
        
        else
        {
            ContigData=readBedGraphFileStrand(BedGraphFile, TrainingIntervals, 1, 0);
        }
    }
    
    else if(p->strand==STRAND_NEGATIVE)
    {
        if(BedGraphSplit)
        {
            ContigData=readSplitBedGraphFileStrand(BedGraphFile, BedGraphNegFile, TrainingIntervals, 1, 1);
        }
        
        else
        {
            ContigData=readBedGraphFileStrand(BedGraphFile, TrainingIntervals, 1, 1);
        }
    }*/
    
    if(strand=="." && (p->strand==STRAND_UNSPECIFIED || p->strand==STRAND_BOTH))
    {
        cout<<"Warning: The input histogram file contains both positive and negative reads and"<<endl;
        cout<<"no specific training strand was specified. This configuration is unfortunately unsupported at this time."<<endl;
        cout<<"Proceeding to train on positive strand data only..."<<endl;
        
        strand="+";
        
        //TODO: Add code to split the input file if it has both types of reads.
        //This will allow us to automatically generate multiple training output files to match strand data.
        
        splitinputfile_t inbeds=splitBedgraphFile(BedGraphFile);
        
        cout<<"Training on positive portion of input histogram..."<<endl;
        
        BedGraphFile=inbeds.in1;
        ContigData=readBedGraphFileAllGivenStrand(BedGraphFile, num_proc, strand);
    
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
            cout<<"RO.EXIT set. exiting..."<<endl;
            return 0;
        }
        if (verbose){
            cout<<"done\n";
            cout<<"Begin parameter estimation, Newtons Method: ";
            cout<<flush;

        }
        //=================================================================
        //NEWTONS METHOD
        vector<double> W;
        try
        {
            //W = learn(RO.X, RO.Y, 0, learning_rate);
            W=ParameterizedNewtonsMethod(RO.X, RO.Y, 0, learning_rate, p->maxConvergenceIters, convergence_threshold);
        } catch(int excep)
        {
            cout<<"Error: regression failed. Exiting...";
            return 0;
        }
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
        writeTrainingFile(outFile+".pos", BWO, learning_rate, max_convergence, convergence_threshold, ChIP, commandLine);
        if (verbose){
            cout<<"done\n";
        }
        
        cout<<"Training on negative portion of input histogram..."<<endl;
        
        BedGraphFile=inbeds.in2;
        strand="-";
        ContigData=readBedGraphFileAllGivenStrand(BedGraphFile, num_proc, strand);
    
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
        RO = run_grabTrainingExamples(R, ContigData, ChIP);
        if (RO.EXIT){
            cout<<"RO.EXIT set. exiting..."<<endl;
            return 0;
        }
        if (verbose){
            cout<<"done\n";
            cout<<"Begin parameter estimation, Newtons Method: ";
            cout<<flush;

        }
        //=================================================================
        //NEWTONS METHOD
        W;
        try
        {
            //W = learn(RO.X, RO.Y, 0, learning_rate);
            W=ParameterizedNewtonsMethod(RO.X, RO.Y, 0, learning_rate, p->maxConvergenceIters, convergence_threshold);
        } catch(int excep)
        {
            cout<<"Error: regression failed. Exiting...";
            return 0;
        }
        if (verbose){
            cout<<"done\n";
            cout<<"Parameter estimation, Baum-Welch          : ";
            cout<<flush;
        }
        //=================================================================
        //BAUM WELCH ALG
        BWO 									= runBW(ContigData, W,max_convergence, convergence_threshold,learning_rate, verbose, num_proc, maxSeed, ChIP);
        if (verbose){
            cout<<"done\n";
            cout<<"Writing learned parameters                : ";
            cout<<flush;

        }
        writeTrainingFile(outFile+".neg", BWO, learning_rate, max_convergence, convergence_threshold, ChIP, commandLine);
        if (verbose){
            cout<<"done\n";
        }
        
        //Delete the split files:
        remove(inbeds.in1ch);
        remove(inbeds.in2ch);
    }
    
    else
    {
        ContigData=readBedGraphFileAllGivenStrand(BedGraphFile, num_proc, strand);
        
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
            cout<<"RO.EXIT set. exiting..."<<endl;
            return 0;
        }
        if (verbose){
            cout<<"done\n";
            cout<<"Begin parameter estimation, Newtons Method: ";
            cout<<flush;

        }
        //=================================================================
        //NEWTONS METHOD
        vector<double> W;
        try
        {
            //W = learn(RO.X, RO.Y, 0, learning_rate);
            W=ParameterizedNewtonsMethod(RO.X, RO.Y, 0, learning_rate, p->maxConvergenceIters, convergence_threshold);
        } catch(int excep)
        {
            cout<<"Error: regression failed. Exiting...";
            return 0;
        }
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
}
