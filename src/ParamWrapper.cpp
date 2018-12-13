#include "ParamWrapper.h"

void ParamWrapper::printUsage()
{
    //We're going to print this block the C way to save on time and pain:
    printf("=================================================================================================\n");
    printf("Fast Stitch Reader (FStitch)\n");
    printf("----------------------------\n");
    printf("A machine learning algorithm tool used to determine regions of active transcription in sequencing a\n");
    printf("  data using log-likelihood regression LLR adjusted hidden Markov Model (HMM).\n\n");
    printf("For more usage information, bug reporting, and questions see the FStitch GitHub repository\n");
    printf("         @Dowell-Lab Organization: https://github.com/Dowell-Lab/FStitch\n");
    printf("==================================================================================================\n\n");
    printf("Usage:                          FStitch [train / segment] [arguments]\n\n");
    printf("==================================================================================================\n\n");
    printf("train                           Trains FStitch on user specified regions in BED4 format.\n");
    printf("                                NOTE: You can only train using one strand (sense or anti-sense).\n\n");
    printf("--------------------------------------------------------------------------------------------------\n\n");
    printf("Required arguments:\n\n");
    printf("  -b || --bedgraph <.bedGraph>   This specifies the bedGraph coverage file.\n");
    printf("  -o || --ouput    <.hmminfo>    This specifies output training file to be used in segment.\n");
    printf("                                 **Must end with the `.hmminfo` extension**.\n");
    // The option "both/." is still written into the code -- however this does not make sense and should be probably removed
    printf("  -s || --strand   <+/->         This specifies whether or not training should be\n");
    printf("                                 performed on either the positive strand, the\n");
    printf("                                 negative strand, or both. The default is 'both'\n");
    printf("  -t || -train     <BED4>        This specifies a training input file containing information for\n");
    printf("                                 BOTH on and off annotations. The input should be a BED3 file\n");
    printf("                                 (chr, start, stop) with an on/off label. For example:\n\n");
    printf("                                 chr [tab] start [tab] stop [tab] 1 / 0 (on / off).\n");
    printf("                                 chr1    2000    7500    1\n");
    printf("                                 chr1    7501    11100   0\n");
    printf("                                 **IMPORTANT: There can be no overlapping training regions.**\n\n");
    printf("Optional arguments:\n\n");
    printf("  -n || --threads <integer>      This specifies the number of threads/processors to run on.\n");
    printf("                                 Default = 1.\n\n");
    /* rp and rn function, but can very easily lead to a numnber of errors if the user inputs a bedGraph with both
     * positive and negative strand information. -r works fine and the strand the user wants to train on based on .bed file
     * should then specified using the argument --strand. NOTE: FStitch cannot train on both strands!!
     */
    //printf(" -rp <pos read BED4>  If -r is not specified, then this specifies the file\n");
    //printf("                      containing the histogram of positive reads as a BED4 file.\n");
    //printf(" -rn <neg read BED4>  If -r is not specified, then this specifies the file\n");
    //printf("                      containing the histogram of negative reads as a BED4 file.\n");
    //printf("                      This parameter must be used in conjunction with -rp\n");
    //printf("                      regardless of the strand used for training.\n");
    /* This flag is somewhat useful (emphasis on somewhat) on our end for trouble-shooting but should not be a documented option
     * until bugs are cleaned up
     */
    //printf(" -v                   This enables verbose logging and output.\n");
    /* This flag appears to still function but leads to many errors -- need to determine how its functioning and either remove
     * or edit it
     */
    //printf(" -chip                This toggles whether or not the input data was\n");
    //printf("                      generated with ChIP-seq or any other single-stranded\n");
    //printf("                      genomic data.\n");
    //TODO: rewrite the following:
    /* Changing this does not appear to have a signficant effect on performance of the model. It's
     * probably best to not have the user specify these
     */
    //printf(" -cm <default=100>    This sets the maximum number of iterations for learning.\n");
    //printf(" -ct <default=0.001>  This sets the convergence threshold.\n");
    //CHANGED FROM -al
    //printf(" -lr <default=0.4>    This sets the learning rate.\n");
    //printf(" -reg <default=1>     This sets some kind of regularization parameter.\n");
    //printf(" -ms <default=20>     This sets the maximum seed value.\n");
    //End of region marked for rewrite.
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
    //printf("You must specify a training file. If --on or --off is specified, then\n");
    //printf("both --on and --off must be specified. -t can be used instead for training\n");
    //printf("files in which both on and off regions are present.\n");
    //printf(" --on <file>       This specifies the training data marking all \"on\"\n");
    //printf("                   regions. It is a bed3 file formatted as follows:\n");
    //printf("                   chromosome[tab]start[tab]stop\n");
    //printf(" --off     <file>  This specifies the training data marking all \"off\"\n");
    //printf("                   regions. It follows the same format as --on.\n");
    printf("==================================================================================================\n\n");
    printf("segment                         Annotates regions of predicted active transcription in BED9 format\n");
    printf("                                using LLR and HMM parameters generated from the train module.\n\n");
    printf("--------------------------------------------------------------------------------------------------\n\n");
    printf("\n");
    printf("Required arguments:\n");
    printf("  -b || --bedgraph <.bedGraph>   This specifies the bedGraph coverage file.\n");
    // The option "both/." is still written into the code -- however this does not make sense and should be probably removed
    printf("  -s || --strand   <+/->         This specifies whether to segment based on information\n");
    printf("                                 in the positive strand, the negative strand, or both.\n");
    printf("                                 This parameter should match what was used in training.\n");
    printf("                                 FStitch will attempt to automatically determine this\n");
    printf("                                 based on the input file provided, but this determination\n");
    printf("                                 may not always be accurate or desirable.\n");
    printf("  -o || --output   <BED9>        This specifies the .bed output annotation file.\n");
    printf("  -p || --params   <.hmminfo>    This specifies the HMM parmaters generated using the \n");
    printf("                                 'train' module.\n\n");
    printf("Optional arguments:\n");
    // This does not currently work -- doesn't even appear to be in anywhere -- probably an easy fix but I don't know how to yet in C world
    //printf("  -r || --report   <on/off/both> This specifies whether annotations generated by FStitch\n");
    //printf("                                 will only report \"on\" regions, \"off\" regions, or both.\n");
    //printf("                                 The default value is \"on\".\n");
    printf("  -n || --threads <integer>     This specifies the number of processors to run on.\n");
    printf("                                 Default = 1\n\n");
    printf("==================================================================================================\n\n");
    printf("  -h || --help                  Prints the help message\n\n");
    // Here again... these arguments don't really make sense in conjunction with the strand argument
    //printf("                        of all reads as a bed4 file.\n");
    //printf("  -rp <pos read BED4>   If -r is not specified, then this specifies the file\n");
    //printf("                        containing the histogram of positive reads as a BED4 file.\n");
    //printf("  -rn <neg read BED4>   If -r is not specified, then this specifies the file\n");
    //printf("                        containing the histogram of negative reads as a BED4 file.\n");
    //printf("                        This parameter must be used in conjunction with -rp\n");
    //printf("                        regardless of the strand used for training.\n");

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
    cout<<"1.1"<<this->version<<endl;
}

//This just sets default parameters:
ParamWrapper::ParamWrapper()
{
    this->commandLine="UNIT TEST FILE";
    this->strand=STRAND_UNSPECIFIED;//STRAND_POSITIVE;
    this->report=REPORT_ON;
    this->exit=false;
    this->version=false;
    this->train=false;
    this->segment=false;
    this->eRNA=false;
    this->verbose=true;
    this->chip=false;
    this->numProcs=1;
    //Set learning parameters:
    this->maxConvergenceIters=100;
    this->learningRate=0.4;
    this->convergenceThreshold=0.001;
    this->regularization=1;
    this->maxSeed=20;
    //Set all split parameters, etc, to blank strings:
    this->specialFileSplit=false;
    this->specialFileName="";
    this->readFileName="";
    this->outFileName="";
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
    this->strand=STRAND_UNSPECIFIED;
    //this->strand=STRAND_BOTH;
    this->report=REPORT_ON;
    this->exit=false;
    this->version=false;
    this->train=false;
    this->segment=false;
    this->eRNA=false;
    this->verbose=false;
    this->chip=false;
    this->numProcs=1;
    
    this->readFileSplit=false;
    this->secondReadFileName="";
    
    // Set some of the additional parameters used to alter how the program learns:
    this->maxConvergenceIters=100;
    this->learningRate=0.4;
    this->convergenceThreshold=0.001;
    this->regularization=1;
    //TODO: find the actual default value for this parameter.
    this->maxSeed=20;
    
    
    if(argc==1)
    {
        printf("ERROR: arguments expected after command specification.\n");
        printf("\n");
        this->exit=true;
        this->printUsage();
        
        return;
    }
    
    if(argc==2)
    {
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

    else
    {
        printf("Error: command not specified. Please specify either 'train' \nor 'segment' immediately after the FStitch command.\n");
        this->printUsage();
        this->exit=true;

        return;
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
            
            else if(!strcmp(prevCmd, "--version"))
            {
                this->version=true;
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
        if(it->first=="-b" || it->first=="--bedgraph")
        {
            this->readFileName=it->second;
            
            if(!fileExists(this->readFileName))
            {
                printf("Error: coverage file specified by -b does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="-rp")
        {
            this->readFileSplit=true;
            this->readFileName=it->second;
            
            //Determine if the specified file exists:
            if(!fileExists(this->readFileName))
            {
                printf("Error: positive coverage file specified by -rp does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="-rn")
        {
            this->readFileSplit=true;
            this->secondReadFileName=it->second;
            
            if(!fileExists(this->secondReadFileName))
            {
                printf("Error: negative coverage file specified by -rn does not exist.\n");
                this->exit=true;
            }
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
        
        else if(it->first=="-o" || it->first=="--output")
        {
            this->outFileName=it->second;
        }
        
        else if(it->first=="-n" || it->first=="threads")
        {
            this->numProcs=atoi(it->second.c_str());
        }
        
        else if(it->first=="--report" || it->first=="-r")
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
        
        else if(it->first=="--strand" || it->first=="-s")
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
        
        else if(it->first=="-t" || it->first=="--train")
        {
            this->specialFileName=it->second;
            // This should take care of the special case in which this value is specified last.
            // It should take priority over individual splits.
            this->specialFileSplit=false;
            
            if(!fileExists(this->specialFileName))
            {
                printf("Error: training input file specified with -t || --train does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="--off")
        {
            this->secondSpecialFileName=it->second;
            this->specialFileSplit=true;
            
            if(!fileExists(this->specialFileName))
            {
                printf("Error: bed file containing off training regions specified with --off does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="--on")
        {
            this->specialFileName=it->second;
            this->specialFileSplit=true;
            
            if(!fileExists(this->specialFileName))
            {
                printf("Error: bed file containing on training regions specified with --on does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="--negtrainingfile")
        {
            this->secondSpecialFileName=it->second;
            this->specialFileSplit=true;
            
            if(!fileExists(this->secondSpecialFileName))
            {
                printf("Error: negative training examples file specified with --negtrainingfile does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="--postrainingfile")
        {
            this->specialFileName=it->second;
            this->specialFileSplit=true;
            
            if(!fileExists(this->specialFileName))
            {
                printf("Error: positive training examples file specified with --postrainingfile does not exist.\n");
                this->exit=true;
            }
        }
        
        else if(it->first=="-p" || it->first=="--params")
        {
            this->specialFileName=it->second;
            this->specialFileSplit=false;
            
            if(!fileExists(this->specialFileName))
            {
                printf("Error: HMM parameters file specified with -p || --params does not exist.\n");
                this->exit=true;
            }
        }
    }
    
    //Set the internal parameters map so that we can kludge our way to better parameter support!
    this->params=paramMap;
    //That should probably do it.
    
    printf("Strand specified (1 = Positive (+), 2 = Negative (-), 3 = both, 4 = unspecified): %d\n", this->strand);
    //If, after all of this, we have an invalid mode:
}
