#include "unittests.h"

using namespace std;

splitoutput_t *genOnOff(string infile)
{
    splitoutput_t *out;
    string line;
    //This nasty mix of C and C++ file I/O brought to you by mkstemp's return value.
    ifstream infile(infile);
    FILE *aFile;
    FILE *bFile;
    
    out=(splitoutput_t *) malloc(sizeof(splitoutput_t));
    
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
    fclose(aFile);
    fclose(bFile);
    return out;
}

splitoutput_t *genPosNeg(string infile)
{
    splitoutput_t *out;
    string line;
    //This nasty mix of C and C++ file I/O brought to you by mkstemp's return value.
    ifstream infile(infile);
    FILE *aFile;
    FILE *bFile;
    
    out=(splitoutput_t *) malloc(sizeof(splitoutput_t));
    
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
        
        if(stoi(toks[3])<0)
        {
            fprintf(bFile, "%s\t%s\t%s\n", toks[0].c_str(), toks[1].c_str(), toks[2].c_str());
        }
        
        else if(stoi(toks[3])>=0)
        {
            fprintf(aFile, "%s\t%s\t%s\n", toks[0].c_str(), toks[1].c_str(), toks[2].c_str());
        }
        
        else
        {
            printf("Unknown label: %s\n", toks[3].c_str());
        }
    }
    
    //These calls also close the corresponding descriptors. This is acceptable behavior.
    fclose(aFile);
    fclose(bFile);
    return out;
}

void cleanupSplitOutput(splitoutput_t *o)
{
    remove(o->aCh);
    remove(o->bCh);
}

void printUnitUsage(char *unitName)
{
    printf("Usage: %s <training bed4 file> <bed4 histogram file>\n");
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
    printf("Generating reference training output...\n");
    w->outFileName="tout_ref.out";
    w->readFileName=bedfile;
    w->specialFileName=trainfile;
    
    main_train(w);
    
    printf("Training with split label file...\n");
    //Now set this up for on/off split training examples:
    w->outFileName="tout_onoffsplit.out";
    w->specialFileName=splitTraining.a;
    w->secondSpecialFileName=splitTraining.b;
    w->specialFileSplit=true;
    
    main_train(w);
    
    //Now just use the set of positive histograms:
    printf("Training with only positive inputs...\n");
    w->outFileName="tout_posonly.out";
    w->specialFileName=trainfile;
    w->specialFileSplit=false;
    //This should be only the positive values:
    w->readFileName=splitBed.a;
    
    main_train(w);
    
    //Now perform segmentation tasks:
    printf("Segmenting with reference inputs...\n");
    
    printf("Segmenting with inputs generated with split label file...\n");
    printf("Segmenting with inputs generated with positive data only...\n");
    
    cleanupSplitOutput(splitTraining);
    cleanupSplitOutput(splitBed);
    free(splitTraining);
    free(splitBed);
    delete w;
    return EXIT_SUCCESS;
}
