#include "main_segment.h"
#include "read.h"
#include <omp.h>
#include <map>
#include "interval_tree.h"
#include "write.h"
#include "viterbi.h"
#include "validate.h"
#include "split.h"

#include <string>
using namespace std;

int run_main_segment_pwrapper(ParamWrapper *p)
{
    //=================================================================
    // General Parameters
    string BedGraphFile 		= p->readFileName;
    string SecondBedGraphFile           = p->secondReadFileName;
    string TrainingOutFile 		= p->specialFileName;
    string outFile 				= p->outFileName;
    //Todo: add functionality for this!
    //string refFile=p->params["-r"];
    string refFile;
    //string refFile 				= PT.params["-r"];
    //This will be handled differently:
    string strand;
    int determinedStrand;

    //==================================================================
    // Other Parameters
    string commandLine=p->commandLine;

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

    int num_proc = p->numProcs;

    //This should be removed:
    /*
    if (strand.empty()){
        if (isPos(BedGraphFile)){
            strand = "+";
        }else if(isNeg(BedGraphFile)){
            strand = "-";
        }else{
            strand = ".";
        }
    }*/
    //
    bool verbose=p->verbose;
    //=================================================================
    //Read in FStich Training out file
    if (verbose){
        cout<<"Reading in FStitch Training Out File: ";
    }
    RTOF RTOF_params 					= readTrainingOutFile(TrainingOutFile);

    if (RTOF_params.EXIT){
        cout<<"RTOF_params.EXIT set. exiting..."<<endl;
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
    
    splitinputfile_t inbeds;
    
    //Were we passed individual bedgraph nuggets?
    if(p->readFileSplit)
    {
        if(verbose)
        {
            cout<<"ParamWrapper reports that input reads file was split."<<endl;
        }
        
        //Determine which value is which:
        if(BedGraphFile=="" && SecondBedGraphFile!="")
        {
            if(strand!="-")
            {
                cerr<<"Error: Mismatch between provided histogram bedgraph and training parameters."<<endl;
                return 0;
            }
            
            else
            {
                inbeds.wasSplit=false;
                inbeds.in1=SecondBedGraphFile;
            }
        }
        
        else if(BedGraphFile!="" && SecondBedGraphFile!="")
        {
            inbeds.wasSplit=true;
            inbeds.in1=BedGraphFile;
            inbeds.in2=SecondBedGraphFile;
            cout<<"in1: "<<inbeds.in1<<" in2: "<<inbeds.in2<<endl;
        }
        
        else if(BedGraphFile!="" && SecondBedGraphFile=="")
        {
            //If this doesn't match the selected strand, then we have a problem:
            if(strand!="+")
            {
                cerr<<"Error: Mismatch between provided histogram bedgraph and training parameters."<<endl;
                return 0;
            }
            
            else
            {
                inbeds.wasSplit=false;
                inbeds.in1=BedGraphFile;
            }
        }
        
        else
        {
            cerr<<"Error: No viable histogram bedgraph file specified. Exiting..."<<endl;
            return 0;
        }
        
        if(verbose)
        {
            cout<<"Passed positive split file name: "<<inbeds.in1<<endl;
            cout<<"Passed negative split file name: "<<inbeds.in2<<endl;
        }
    }
    
    else
    {
        if(verbose)
        {
            cout<<"Splitting input bed file."<<endl;
        }
        
        //It may be the case that the input bedgraph was improperly specified.
        //isPos and isNeg should've made this determination earlier:
        if(strand=="+" || strand=="-")
        {
            inbeds.wasSplit=false;
            
            if(verbose)
            {
                cout<<"Histogram file contains only histograms for one strand."<<endl;
            }
        }
        
        else
        {
            inbeds=splitBedgraphFile(BedGraphFile);
            
            if(verbose)
            {
                cout<<"Positive split file name: "<<inbeds.in1<<endl;
                cout<<"Negative split file name: "<<inbeds.in2<<endl;
            }
        }
    }

    cout<<"Input bed file split? "<<inbeds.wasSplit<<endl;

    //In this case, we need to segment twice:
    if(inbeds.wasSplit)
    {
        vector<string> outFilePathToks=splitter(outFile, "/");
        string outFilePathPrefix="";
        for(int i=0;i<outFilePathToks.size()-1;i++)
        {
            outFilePathPrefix+=outFilePathToks[i]+"/";
        }
        
        vector<string> outFileToks=splitter(outFilePathToks[outFilePathToks.size()-1], ".");

        map<string,contig *> ContigData = readBedGraphFileAllGivenStrand(inbeds.in1, num_proc, "+", verbose);
        if (ContigData.empty()){
            cout<<"ContigData is empty. exiting..."<<endl;
            return 0;
        }

        if (verbose){
            cout<<"done"<<endl;
            cout<<flush;
        }
        //=================================================================
        //Run Viterbi
        if (verbose){
            cout<<"Running Viterbi on positive examples        : ";
            cout<<flush;
        }
        map<string, state*> results = runViterbi(ContigData, RTOF_params.W, RTOF_params.A,num_proc, RTOF_params.ChIP);
        if (verbose){
            cout<<"done"<<endl;
            cout<<flush;
        }
        //=================================================================
        //Write Viterbi Paths
        if (verbose){
            cout<<"Writing positive to IGV                      : ";
            cout<<flush;
        }
        
        writeViterbiPaths(outFilePathPrefix+outFileToks[0]+".pos.bed", results, refFile, strand, RTOF_params.commandLine, commandLine);
        if (verbose){
            cout<<"done"<<endl;
        }
        
        map<string, contig *>::iterator cti;
        int ctgChainSize=0;
        
        for(cti=ContigData.begin(); cti!=ContigData.end();cti++)
        {
            if(cti->second)
            {
                ctgChainSize+=altDelContigChain(cti->second, 1);
                //delete(cti->second);
                cti->second=NULL;
            }
        }
        
        ContigData.clear();
        
        printf("Deleted a %d contigs totaling %lu bytes.\n", ctgChainSize, ctgChainSize*sizeof(contig));

        ContigData 	= readBedGraphFileAllGivenStrand(inbeds.in2, num_proc, "-", verbose);
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
            cout<<"Running Viterbi on negative examples         : ";
            cout<<flush;
        }

        results 	= runViterbi(ContigData, RTOF_params.W, RTOF_params.A,num_proc, RTOF_params.ChIP);
        if (verbose){
            cout<<"done"<<endl;
            cout<<flush;
        }
        //=================================================================
        //Write Viterbi Paths
        if (verbose){
            cout<<"Writing negative to IGV                      : ";
            cout<<flush;
        }
        writeViterbiPaths(outFilePathPrefix+outFileToks[0]+".neg.bed", results, refFile, strand, RTOF_params.commandLine, commandLine);
        if (verbose){
            cout<<"done"<<endl;
        }
        
        ctgChainSize=0;
        
        for(cti=ContigData.begin(); cti!=ContigData.end();cti++)
        {
            if(cti->second)
            {
                ctgChainSize+=altDelContigChain(cti->second, 1);
                //delete(cti->second);
                cti->second=NULL;
            }
        }
        
        ContigData.clear();
        
        printf("Deleted a %d contigs totaling %lu bytes.\n", ctgChainSize, ctgChainSize*sizeof(contig));

        //Only delete the files if they were either generated here or we don't want debugging information.
        if(!verbose && inbeds.splitGenerated)
        {
            remove(inbeds.in1ch);
            remove(inbeds.in2ch);
        }
        
        else
        {
            cout<<"All temporary files mentioned previously have been preserved. Delete them manually if necessary."<<endl;
        }
    }

    //Otherwise, we just need to segment once:
    else
    {
        map<string,contig *> ContigData 	= readBedGraphFileAllGivenStrand(BedGraphFile, num_proc, strand, verbose);//readBedGraphFileAll(BedGraphFile,num_proc);
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

        //Now delete the input files:
        if(!verbose && inbeds.splitGenerated)
        {
            remove(inbeds.in1ch);
        }
        
        else
        {
            cout<<"All temporary files mentioned previously have been preserved. Delete them manually if necessary."<<endl;
        }
        
        map<string, contig *>::iterator cti;
        int ctgChainSize=0;
        
        for(cti=ContigData.begin(); cti!=ContigData.end();cti++)
        {
            if(cti->second)
            {
                ctgChainSize+=altDelContigChain(cti->second, 1);
                //delete(cti->second);
                cti->second=NULL;
            }
        }
        
        ContigData.clear();
        
        printf("Deleted a %d contigs totaling %lu bytes.\n", ctgChainSize, ctgChainSize*sizeof(contig));
    }

    return 1;
}
