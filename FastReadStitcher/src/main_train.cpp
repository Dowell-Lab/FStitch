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
bool EXIT(string BedGraphFile, string TrainingFile, string outFile, string np, string mc, string lr, string ct){
	bool BGF 	= isFile(BedGraphFile);
	bool TOF 	= isFile(TrainingFile);
	bool NP 	= isNum(np);
	bool MC 	= isNum(mc);
	bool LR 	= isNum(lr);
	bool MS 	= isNum(ct);
	if (not BGF){
		cout<<"FILE (-i): "<<BedGraphFile<<", does not exist"<<endl;
		return 1;
	}
	if (not TOF){
		cout<<"Training File (-j): "<< TrainingFile<<", does not exist"<<endl;
		return 1;
	}
	if (not NP){
		cout<<"Number of Processors (-np): "<<np<<" is not an integer value"<<endl;
		return 1;
	}
	if (not MC){
		cout<<"Maximum number of learning iterations (-cm): "<<mc<<" is not a number"<<endl;
		return 1;
	}
	if (not LR){
		cout<<"Learning rate (-lr): "<<lr<<" is not a number"<<endl;
		return 1;
	}
	if (not MS){
		cout<<"Convergence threshold (-ct): "<<ct<<" is not a number"<<endl;
		return 1;
	}

	return 0;

}

int run_main_train(paramsTrain PT){
	bool exit_bool 								= EXIT(PT.params["-i"], PT.params["-j"], PT.params["-o"], PT.params["-np"],
		PT.params["-cm"], PT.params["-lr"], PT.params["-ct"]);
	if (exit_bool){
		cout<<"exiting..."<<endl;
		return 0;
	}

	//=================================================================
	// General Parameters

	string BedGraphFile 						= PT.params["-i"];
	string TrainingFile 						= PT.params["-j"];
	string outFile 								= PT.params["-o"];
	int num_proc 								= stoi(PT.params["-np"]);
	bool verbose 								= not PT.params["-v"].empty();
	bool ChIP 									= not PT.params["-chip"].empty();
	//=================================================================
	// parameters specific to Baum Welch and Newtons Method
	double max_convergence 						= stof(PT.params["-cm"]);
	double convergence_threshold 				= stof(PT.params["-ct"]);
	double learning_rate 						= stof(PT.params["-al"]);
	int maxSeed 								= stoi(PT.params["-ms"]);
	
	//=================================================================

	if (verbose){
		cout<<"reading training file                     : ";
	}
	//=================================================================
	//READ TRAINING FILE
	readTrainingFileReturn TrainReturn 	= readTrainingFile(TrainingFile);
	if (TrainReturn.EXIT){
		cout<<"exiting..."<<endl;
		return 0;
	}
	map<string, interval *> TrainingIntervals 	= TrainReturn.result;

	if (verbose){
		cout<<"done"<<endl;
		cout<<"making interval tree                      : ";
	}
	//=================================================================
	//INTERVAL TREE FROM TRAINING FILE
	map<string,T> R 							= makeIntervalTree(TrainingIntervals);
	if (verbose){
		cout<<"done"<<endl;
		cout<<"reading bed graph file                    : ";
		cout<<flush;
	}
	//=================================================================
	//READ BEDGRAPH FILE
	map<string,contig *> ContigData 			= readBedGraphFile(BedGraphFile,TrainingIntervals,1);
	if (ContigData.empty()){
		cout<<"exiting..."<<endl;
		return 0;
	}
	if (verbose){
		cout<<"done\n";
		cout<<"grabbing training data from bedgraph file : ";
		cout<<flush;
	
	}
	//=================================================================
	//GET DATA FROM TRAINING INTERVALS 
	run_out RO  								= run_grabTrainingExamples(R, ContigData, ChIP);
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
	writeTrainingFile(outFile, BWO, learning_rate, max_convergence, convergence_threshold, ChIP);
	if (verbose){
		cout<<"done\n";
    }
}	

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
    writeTrainingFile(outFile, BWO, learning_rate, max_convergence, convergence_threshold, ChIP);
    if (verbose){
        cout<<"done\n";
    }
}
