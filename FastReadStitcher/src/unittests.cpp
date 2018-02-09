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
    remove(o->aCh);
    remove(o->bCh);
}

void printUnitUsage(char *unitName)
{
    printf("Usage: %s <training bed4 file> <bed4 histogram file>\n", unitName);
    printf("Runs unit tests for FStitch given input training and histogram files.\n");
}

/* This value represents how far "off" a given value can be from */
#define ACCEPTABLE_MARGIN 0.01

bool checkWeightsConsistency(string goldStandard, string checkFile, double errorMargin)
{
    vector<string> goldStandardLines;
    vector<string> checkFileLines;
    vector<string> gsToks;
    vector<string> cfToks;
    int goldStandardStart;
    int checkStart;
    int i, j;
    int nRemaining;
    double g, c;
    
    goldStandardLines=readFileLines(goldStandard);
    checkFileLines=readFileLines(checkFile);
    
    goldStandardStart=0;
    checkStart=0;
    
    while(goldStandardLines[goldStandardStart][0]=='#')
    {
        goldStandardStart++;
    }
    
    while(checkFileLines[checkStart][0]=='#')
    {
        checkStart++;
    }
    
    //Now that we have appropriate starting positions, let's attempt to compare our outputs.
    //This first comparison should be whether or not the training process converged:
    gsToks=splitter(goldStandardLines[goldStandardStart], ":");
    cfToks=splitter(checkFileLines[checkStart], ":");
    
    if(gsToks[1]!=cfToks[1])
    {
        return false;
    }
    goldStandardStart++;
    checkStart++;
    
    //The second comparison is of the final log likelihood:
    g=atof(splitter(goldStandardLines[goldStandardStart], ":")[1].c_str());
    c=atof(splitter(checkFileLines[checkStart], ":")[1].c_str());
    
    if(fabs(g-c)>errorMargin)
    {
        return false;
    }
    
    goldStandardStart++;
    checkStart++;
    
    nRemaining=goldStandardLines.size()-goldStandardStart;
    
    //Now we get to the fun part. Each of the next few lines contains a comma separated list of values.
    for(i=0;i<nRemaining;i++)
    {
        //This is an awful implementation;
        gsToks=splitter(splitter(goldStandardLines[goldStandardStart+i], ":")[1], ",");
        cfToks=splitter(splitter(checkFileLines[checkStart+i], ":")[1], ","); 
        
        if(gsToks.size()!=cfToks.size())
        {
            printf("Line size mismatch for gold standard %d and check %d\n", goldStandardStart+i, checkStart+i);
            return false;
        }
        
        for(j=0;j<gsToks.size();j++)
        {
            g=atof(gsToks[j].c_str());
            c=atof(cfToks[j].c_str());
            
            if(fabs(g-c)>errorMargin)
            {
                return false;
            }
        }
    }
    
    return true;
}

bool checkBedsConsistency(string goldStandard, string checkFile)
{
    vector<string> goldStandardLines;
    vector<string> checkFileLines;
    int i;
    
    goldStandardLines=readFileLines(goldStandard);
    checkFileLines=readFileLines(checkFile);
    
    if(goldStandardLines.size()!=checkFileLines.size() && false)
    {
        printf("Inconsistent number of lines: %lu vs %lu.\n", goldStandardLines.size(), checkFileLines.size());
        return false;
    }
    
    for(i=1;i<goldStandardLines.size();i++)
    {
        if(goldStandardLines[i]!=checkFileLines[i])
        {
            printf("Inconsistency found on line %d.\n", i+1);
            printf("Gold standard line: %s\n", goldStandardLines[i].c_str());
            printf("Check file line: %s\n", checkFileLines[i].c_str());
            return false;
        }
    }
    
    return true;
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
    w->verbose=true;
    
    //Generate reference data with default parameters (relatively speaking):
    
    printf("\nGenerating reference training output...\n");
    w->outFileName="tout_ref.out";
    w->readFileName=bedfile;
    w->specialFileName=trainfile;
    
    w->dumpValues();
    
    run_main_train_pwrapper(w);
    
    
    printf("\nTraining with split label file...\n");
    //Now set this up for on/off split training examples:
    w->outFileName="tout_onoffsplit.out";
    w->specialFileName=splitTraining->a;
    w->secondSpecialFileName=splitTraining->b;
    w->specialFileSplit=true;
    
    run_main_train_pwrapper(w);
    
    //Now use a split input histogram:
    printf("\nTraining with split input histogram bedgraph...\n");
    w->outFileName="tout_splitbed.out";
    w->specialFileName=trainfile;
    w->specialFileSplit=false;
    //This should be only the positive values:
    w->readFileSplit=true;
    w->readFileName=splitBed->a;
    w->secondReadFileName=splitBed->b;
    run_main_train_pwrapper(w);
    
    //Now split both types of inputs:
    w->outFileName="tout_bothsplit.out";
    w->specialFileName=splitTraining->a;
    w->secondSpecialFileName=splitTraining->b;
    w->specialFileSplit=true;
    run_main_train_pwrapper(w);
    
    printf("\nDetermining if all training output files are consistent...\n");
    printf("tout_ref.out vs tout_onoffsplit.out.........");
    checkWeightsConsistency("tout_ref.out", "tout_onoffsplit.out", ACCEPTABLE_MARGIN) ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.out vs tout_splitbed.out...........");
    checkWeightsConsistency("tout_ref.out", "tout_splitbed.out", ACCEPTABLE_MARGIN) ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.out vs tout_bothsplit.out..........");
    checkWeightsConsistency("tout_ref.out", "tout_bothsplit.out", ACCEPTABLE_MARGIN) ? printf("PASS\n") : printf("FAIL\n");
    
    //Now perform segmentation tasks:
    
    printf("\nSegmenting with reference inputs...\n");
    w->outFileName="tout_ref.bed";
    //w->strand=STRAND_POSITIVE;
    w->strand=STRAND_UNSPECIFIED;
    w->readFileSplit=false;
    w->readFileName=bedfile;
    //This corresponds to the -w parameter:
    w->specialFileName="tout_ref.out";
    w->specialFileSplit=false;
    run_main_segment_pwrapper(w);
    
    printf("\nSegmenting with inputs generated with split label file...\n");
    w->outFileName="tout_onoffsplit.bed";
    w->specialFileName="tout_onoffsplit.out";
    run_main_segment_pwrapper(w);
    
    printf("\nSegmenting with inputs generated with split input histogram bedgraph...\n");
    w->outFileName="tout_splitbed.bed";
    w->specialFileName="tout_splitbed.out";
    run_main_segment_pwrapper(w);
    
    //printf("\nSegmenting with inputs generated with only positive points on input histogram...\n");
    //w->outFileName="tout_posbed.bed";
    //w->specialFileName="tout_posbed.bed";
    //run_main_segment_pwrapper(w);
    
    printf("\nSegmenting with reference inputs given split histogram.\n");
    w->outFileName="tout_splithist.bed";
    w->specialFileName="tout_ref.out";
    w->readFileSplit=true;
    w->readFileName=splitBed->a;
    w->secondReadFileName=splitBed->b;
    run_main_segment_pwrapper(w);
    
    //Now that we have output files, we need to ensure that they're all consistent:
    printf("\nDetermining if all segmentation output files are consistent given the same training inputs...\n");
    printf("tout_ref.pos.bed vs tout_onoffsplit.pos.bed.....");
    checkBedsConsistency("tout_ref.pos.bed", "tout_onoffsplit.pos.bed") ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.neg.bed vs tout_onoffsplit.neg.bed.....");
    checkBedsConsistency("tout_ref.neg.bed", "tout_onoffsplit.neg.bed") ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.pos.bed vs tout_splitbed.pos.bed.......");
    checkBedsConsistency("tout_ref.pos.bed", "tout_splitbed.pos.bed") ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.neg.bed vs tout_splitbed.neg.bed........");
    checkBedsConsistency("tout_ref.neg.bed", "tout_splitbed.neg.bed") ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.pos.bed vs tout_splithist.pos.bed.......");
    checkBedsConsistency("tout_ref.pos.bed", "tout_splithist.pos.bed") ? printf("PASS\n") : printf("FAIL\n");
    printf("tout_ref.neg.bed vs tout_splithist.neg.bed.......");
    checkBedsConsistency("tout_ref.neg.bed", "tout_splithist.neg.bed") ? printf("PASS\n") : printf("FAIL\n");
    
    cleanupSplitOutput(splitTraining);
    cleanupSplitOutput(splitBed);
    free(splitTraining);
    free(splitBed);
    delete w;
    return EXIT_SUCCESS;
}
