#include "unittests.h"

using namespace std;

splitoutput_t *genOnOff(string infilename)
{
    splitoutput_t *out;
    string line;
    //This nasty mix of C and C++ file I/O brought to you by mkstemp's return value.
    ifstream infile(infilename);
    FILE *aFile;
    FILE *bFile;
    
    out=(splitoutput_t *) malloc(sizeof(splitoutput_t));
    //This is simultaneously for safety and to shut valgrind up about uninitialized values:
    memset((void*) out, 0, sizeof(splitoutput_t));
    
    //Open applicable file descriptors:
    strcpy(out->aCh, "/tmp/tmpXXXXXX");
    out->aFd=mkstemp(out->aCh);
    if(out->aFd<0)
    {
        printf("Error: could not open temporary file %s\n", out->aCh);
    }
    out->a=out->aCh;
    aFile=fdopen(out->aFd, "w");
    
    strcpy(out->bCh, "/tmp/tmpXXXXXX");
    out->bFd=mkstemp(out->bCh);
    if(out->bFd<0)
    {
        printf("Error: could not open temporary file %s\n", out->bCh);
    }
    out->b=out->bCh;
    bFile=fdopen(out->bFd, "w");
    
    while(getline(infile, line))
    {
        vector<string> toks=splitter(line, "\t");
        
        if(toks[3]=="0")
        {
            fprintf(bFile, "%s\t%s\t%s\n", toks[0].c_str(), toks[1].c_str(), toks[2].c_str());
        }
        
        else if(toks[3]=="1")
        {
            fprintf(aFile, "%s\t%s\t%s\n", toks[0].c_str(), toks[1].c_str(), toks[2].c_str());
        }
        
        else
        {
            printf("Unknown label: %s\n", toks[3].c_str());
        }
    }
    
    //These calls also close the corresponding descriptors. This is acceptable behavior.
    //fclose(aFile);
    //fclose(bFile);
    out->aFile=aFile;
    out->bFile=bFile;
    
    fflush(aFile);
    fflush(bFile);
    return out;
}

splitoutput_t *genPosNeg(string infilename)
{
    splitoutput_t *out;
    string line;
    //This nasty mix of C and C++ file I/O brought to you by mkstemp's return value.
    ifstream infile(infilename);
    FILE *aFile;
    FILE *bFile;
    
    out=(splitoutput_t *) malloc(sizeof(splitoutput_t));
    memset((void*) out, 0, sizeof(splitoutput_t));
    
    //Open applicable file descriptors:
    strcpy(out->aCh, "/tmp/tmpXXXXXX");
    out->aFd=mkstemp(out->aCh);
    out->a=out->aCh;
    aFile=fdopen(out->aFd, "w");
    
    strcpy(out->bCh, "/tmp/tmpXXXXXX");
    out->bFd=mkstemp(out->bCh);
    out->b=out->bCh;
    bFile=fdopen(out->bFd, "w");
    
    while(getline(infile, line))
    {
        vector<string> toks=splitter(line, "\t");
        
        //We still expect bed4 files regardless. TODO: Add support for split bedgraph histogram data.
        if(stoi(toks[3])<0)
        {
            fprintf(bFile, "%s\t%s\t%s\t%s\n", toks[0].c_str(), toks[1].c_str(), toks[2].c_str(), toks[3].c_str());
        }
        
        else if(stoi(toks[3])>=0)
        {
            fprintf(aFile, "%s\t%s\t%s\t%s\n", toks[0].c_str(), toks[1].c_str(), toks[2].c_str(), toks[3].c_str());
        }
        
        else
        {
            printf("Unknown label: %s\n", toks[3].c_str());
        }
    }
    
    //These calls also close the corresponding descriptors. This is acceptable behavior.
    //fclose(aFile);
    //fclose(bFile);
    out->aFile=aFile;
    out->bFile=bFile;
    
    fflush(aFile);
    fflush(bFile);
    return out;
}

void cleanupSplitOutput(splitoutput_t *o)
{
    fclose(o->aFile);
    fclose(o->bFile);
    //remove(o->aCh);
    //remove(o->bCh);
}

void printUnitUsage(char *unitName)
{
    printf("Usage: %s <training bed4 file> <bed4 histogram file>\n", unitName);
    printf("Runs unit tests for FStitch given input training and histogram files.\n");
}

int main(int argc, char **argv)
{
    ParamWrapper *w;
    splitoutput_t *splitTraining, *splitBed;
    
    /* Our parameters are fixed: */
    if(argc!=3)
    {
        printUnitUsage(argv[0]);
        return EXIT_SUCCESS;
    }
    
    //This hurts as a pure C programmer:
    string trainfile(argv[1]);
    string bedfile(argv[2]);
    
    splitTraining=genOnOff(trainfile);
    splitBed=genPosNeg(bedfile);
    
    w=new ParamWrapper();
    
    //Generate reference data with default parameters (relatively speaking):
    /*
    printf("\nGenerating reference training output...\n");
    w->outFileName="tout_ref.out";
    w->readFileName=bedfile;
    w->specialFileName=trainfile;
    
    w->dumpValues();
    
    run_main_train_pwrapper(w);
    */
    /*
    printf("\nTraining with split label file...\n");
    //Now set this up for on/off split training examples:
    w->outFileName="tout_onoffsplit.out";
    w->specialFileName=splitTraining->a;
    w->secondSpecialFileName=splitTraining->b;
    w->specialFileSplit=true;
    
    run_main_train_pwrapper(w);
    */
    //Now just use the set of positive histograms:
    printf("\nTraining with only positive inputs...\n");
    w->outFileName="tout_posonly.out";
    w->specialFileName=trainfile;
    w->specialFileSplit=false;
    //This should be only the positive values:
    w->readFileSplit=true;
    w->readFileName=splitBed->a;
    w->secondReadFileName=splitBed->b;
    
    
    run_main_train_pwrapper(w);
    
    //Now perform segmentation tasks:
    printf("\nSegmenting with reference inputs...\n");
    
    printf("\nSegmenting with inputs generated with split label file...\n");
    printf("\nSegmenting with inputs generated with positive data only...\n");
    
    cleanupSplitOutput(splitTraining);
    cleanupSplitOutput(splitBed);
    free(splitTraining);
    free(splitBed);
    delete w;
    return EXIT_SUCCESS;
}
