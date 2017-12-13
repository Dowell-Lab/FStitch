#include "ParamWrapper.h"

void ParamWrapper::printUsage()
{
    //We're going to print this block the C way to save on time and pain:
    printf("Usage: FStitch command [arguments]\n");
    printf("Where command is one of the following:\n");
    printf(" train   -- Trains FStitch on a set of training data specified using\n");
    printf("            some combination of arguments.\n");
    printf(" segment -- Generates annotations based on a trained model and input\n");
    printf("            data specified using some combination of arguments.\n");
    printf("\n");
    printf("If 'train' is specified, then [arguments] may be one or more of the\n");
    printf("following:\n");
    printf(" -r <read bedgraph>   This specifies the file containing  the histogram\n");
    printf("                      of all reads.\n");
    printf(" -o <output file>     This specifies the file to store.\n");
    printf(" -np <number>         This specifies the number of processors to run on.\n");
    printf("                      The default value is 8.\n");
    printf(" -v                   This enables verbose logging and output.\n");
    printf(" -chip                This toggles whether or not the input data was\n");
    printf("                      generated with ChIP-seq or any other single-stranded\n");
    printf("                      genomic data.\n");
    //TODO: rewrite the following:
    printf(" -cm <default=100>    This sets the maximum number of iterations for learning.\n");
    printf(" -ct <default=0.001>  This sets the convergence threshold.\n");
    //CHANGED FROM -al
    printf(" -lr <default=0.4>    This sets the learning rate.\n");
    printf(" -reg <default=1>     This sets some kind of regularization parameter.\n");
    printf(" -ms <default=20>     This sets the maximum seed value.\n");
    //End of region marked for rewrite.
    printf(" --strand <+/-/both>  This specifies whether or not training should be\n");
    printf("                      performed on either the positive strand, the\n");
    printf("                      negative strand, or both. The default is 'both'\n");
    //printf("If both strands are selected:\n");
    //printf(" --postrainingfile <file> This specifies the training data for the positive strand.\n");
    //printf("                   This file is a bed3 file with labels. Ie. it uses the\n");
    //printf("                   following format:\n");
    //printf("                   chromosome[tab]start[tab]stop[tab]1 or 0 for on/off.\n");
    //printf(" --negtrainingfile <file> This specifies the training data for the negative strand.\n");
    //printf("                   This file should be formatted the same as the argument to\n");
    //printf("                   posfile.\n");
    //printf(" -a <file>         This specifies an annotations file containing information\n");
    //printf("                   for BOTH positive and negative strands. If this option is\n");
    //printf("                   specified, then neither posfile nor negfile need be used.\n");
    //printf("                   This file is a bed6 file with the following specification:\n");
    //printf("                   chr[tab]start[tab]stop[tab]1 or 0[tab]score[tab]strand\n");
    //printf("\n");
    //printf("If only a single strand was selected:\n");
    printf("You must speciy a training file. If --on or --off is specified, then\n");
    printf("both --on and --off must be specified. -t can be used instead for training\n");
    printf("files in which both on and off regions are present.\n");
    printf(" --on <file>       This specifies the training data marking all \"on\"\n");
    printf("                   regions. It is a bed3 file formatted as follows:\n");
    printf("                   chromosome[tab]start[tab]stop\n");
    printf(" --off     <file>  This specifies the training data marking all \"off\"\n");
    printf("                   regions. It follows the same format as --on.\n");
    printf(" -t <file>         This specifies a training input file containing information\n");
    printf("                   for BOTH on and off annotations. If this option is specified,\n");
    printf("                   then neither onfile nor negfile need be used. The input should\n");
    printf("                   be a bed3 file with labels, ie. it uses the following format:\n");
    printf("                   chr[tab]start[tab]stop[tab]1 or 0 for on/off.\n");
    printf("\n");
    printf("If 'segment' is specified, then [arguments] may be one or more of the\n");
    printf("following:\n");
    printf(" --strand <+/->         This specifies whether to segment based on information\n");
    printf("                        in the positive strand, the negative strand, or both.\n");
    printf("                        This parameter should match what was used in training.\n");
    printf("                        FStitch will attempt to automatically determine this\n");
    printf("                        based on the input file provided, but this determination\n");
    printf("                        may not always be accurate or desirable.\n");
    printf(" --report <on/off/both> This specifies whether annotations generated by FStitch\n");
    printf("                        should only report \"on\" regions, \"off\" regions, or both.\n");
    printf("                        The default value is \"on\".\n");
    printf(" -w <file>              This specifies the input weights generated using the 'train'\n");
    printf("                        command.\n");

    printf(" -r <read bedgraph>     This specifies the file containing  the histogram\n");
    printf("                        of all reads.\n");
    printf(" -o <file>              This specifies the output annotation file.\n");
    printf("");
}

void ParamWrapper::dumpValues()
{
    //We're going to use C++ output here to avoid mixing printf and C++ strings:
    cout<<"Parameter dump:"<<endl;
    cout<<"outFileName: "<<this->outFileName<<endl;
    cout<<"specialFileName: "<<this->specialFileName<<endl;
    cout<<"readFileName: "<<this->readFileName<<endl;
    cout<<"specialFileSplit: "<<this->specialFileSplit<<endl;
    cout<<"secondSpecialFileName: "<<this->secondSpecialFileName<<endl;
    cout<<"strand: "<<this->strand<<endl;
    cout<<"report: "<<this->report<<endl;
    cout<<"train: "<<this->train<<endl;
    cout<<"segment: "<<this->segment<<endl;
    cout<<"exit: "<<this->exit<<endl;
}

ParamWrapper::ParamWrapper(int argc, char **argv)
{
    int i;
    map<string, string> paramMap;
    char *prevCmd;
    map<string, string>::iterator it;
    
    this->commandLine="";
    
    for(i=0;i<argc;i++)
    {
        commandLine+=argv[i];
        commandLine+=" ";
    }
    
    // Set default parameters:
    //Changed to STRAND_POSITIVE to work around bugs in contig generation.
    this->strand=STRAND_POSITIVE; //STRAND_UNSPECIFIED;
    //this->strand=STRAND_BOTH;
    this->report=REPORT_ON;
    this->exit=false;
    this->train=false;
    this->segment=false;
    this->eRNA=false;
    this->verbose=false;
    this->chip=false;
    this->numProcs=8;
    
    // Set some of the additional parameters used to alter how the program learns:
    this->maxConvergenceIters=100;
    this->learningRate=0.4;
    this->convergenceThreshold=0.001;
    this->regularization=1;
    //TODO: find the actual default value for this parameter.
    this->maxSeed=20;
    
    
    if(argc==1)
    {
        this->exit=true;
        this->printUsage();
        
        return;
    }
    
    if(argc==2)
    {
        printf("Error: arguments expected after command specification.\n");
        this->exit=true;
        this->printUsage();
        
        return;
    }
    
    // Determine what the command is:
    if(!strcmp(argv[1], "segment"))
    {
        this->segment=true;
    }
    
    else if(!strcmp(argv[1], "train"))
    {
        this->train=true;
    }
    
    // Segment the inputs into argument, value pairs:
    for(i=2;i<argc;i++)
    {
        if(argv[i][0]=='-' && strlen(argv[i])!=1)
        {
            prevCmd=argv[i];
            //Check to ensure that the argument is supposed to be able to work without additional parameters:
            if(!strcmp(prevCmd, "-v"))
            {
                this->verbose=true;
            }
            
            else if((!strcmp(prevCmd, "-h"))||(!strcmp(prevCmd, "--help")))
            {
                this->printUsage();
                this->exit=true;
                return;
            }
            
            else if(!strcmp(prevCmd, "-chip"))
            {
                this->chip=true;
            }
        }
        
        else
        {
            paramMap.insert(pair<string, string>(prevCmd, argv[i]));
            prevCmd=NULL;
        }
    }
    
    // Now that we have a map of parameters:
    for(it=paramMap.begin();it!=paramMap.end();it++)
    {
        // These comparisons should be acceptable, since the C++ string library
        //overrides the == operator.
        if(it->first=="-r")
        {
            this->readFileName=it->second;
        }
        
        else if(it->first=="-cm")
        {
            this->maxConvergenceIters=atoi(it->second.c_str());
        }
        
        else if(it->first=="-ct")
        {
            this->convergenceThreshold=atof(it->second.c_str());
        }
        
        else if(it->first=="-lr")
        {
            this->learningRate=atof(it->second.c_str());
        }
        
        else if(it->first=="-reg")
        {
            this->regularization=atoi(it->second.c_str());
        }
        
        else if(it->first=="-o")
        {
            this->outFileName=it->second;
        }
        
        else if(it->first=="-np")
        {
            this->numProcs=atoi(it->second.c_str());
        }
        
        else if(it->first=="--report")
        {
            if(it->second=="on")
            {
                this->report=REPORT_ON;
            }
            
            else if(it->second=="off")
            {
                this->report=REPORT_OFF;
            }
            
            else if(it->second=="both")
            {
                this->report=REPORT_BOTH;
            }
            
            else
            {
                printf("Error: %s is not a legal operand for %s.\n", it->second.c_str(), it->first.c_str());
                this->exit=true;
                this->printUsage();
                
                return;
            }
        }
        
        else if(it->first=="--strand")
        {
            if(it->second=="+")
            {
                this->strand=STRAND_POSITIVE;
            }
            
            else if(it->second=="-")
            {
                this->strand=STRAND_NEGATIVE;
            }

            else if(it->second=="both" || it->second==".")
            {
                this->strand=STRAND_BOTH;
            }
            
            else
            {
                printf("Error: %s is not a legal operand for %s.\n", it->second.c_str(), it->first.c_str());
                this->exit=true;
                this->printUsage();
                
                return;
            }
        }
        
        else if(it->first=="-a")
        {
            this->specialFileName=it->second;
            // This should take care of the special case in which this value is specified last.
            // It should take priority over individual splits.
            this->specialFileSplit=false;
        }
        
        else if(it->first=="--off")
        {
            this->secondSpecialFileName=it->second;
        }
        
        else if(it->first=="--on")
        {
            this->specialFileName=it->second;
            this->specialFileSplit=true;
        }
        
        else if(it->first=="--negtrainingfile")
        {
            this->secondSpecialFileName=it->second;
        }
        
        else if(it->first=="--postrainingfile")
        {
            this->specialFileName=it->second;
            this->specialFileSplit=true;
        }
        
        else if(it->first=="-w")
        {
            this->specialFileName=it->second;
            this->specialFileSplit=false;
        }
    }
    
    //Set the internal parameters map so that we can kludge our way to better parameter support!
    this->params=paramMap;
    //That should probably do it.
    
    printf("Strand specified: %d\n", this->strand);
}
