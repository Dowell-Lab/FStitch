#include "read.h"
#include <string>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <map>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "BaumWelch.h"
#include "split.h"
#include "validate.h"
#include <stdexcept>
using namespace std;

splitinputfile_t splitBedgraphFile(string infilename)
{
    ifstream infile(infilename);
    string line;
    vector<string> lineToks;
    splitinputfile_t out;
    ofstream outfile1, outfile2;

    strcpy(out.in1ch, "/tmp/tmpXXXXXX");
    strcpy(out.in2ch, "/tmp/tmpXXXXXX");

    //TODO: Figure out what to do with the file descriptor that this function returns.
    //We just want the filename, not the handle.
    //Should I close it?
    close(mkstemp(out.in1ch));
    out.in1=string(out.in1ch);
    outfile1=ofstream(out.in1);
    
    out.splitGenerated=true;
    out.wasSplit=false;

    while(getline(infile, line))
    {
        lineToks=splitter(line, "\t");

        //If we're negative:
        if(stoi(lineToks[3])<0)
        {
            if(!out.wasSplit)
            {
                out.wasSplit=true;
                close(mkstemp(out.in2ch));
                out.in2=string(out.in2ch);
                outfile2=ofstream(out.in2);
            }

            //We need to write the tokens individually:
            outfile2<<line<<endl; //lineToks[0]<<"\t"<<lineToks[1]<<"\t"<<lineToks[2]<<"\t"<<stoi(lineToks[3])*-1<<endl;
        }

        //If we're positive:
        else
        {
            outfile1<<line<<endl;
        }
    }

    return out;
}

string getChrom(string line){
	const char * tab = "\t";
	for (int i = 0; i < line.size(); i++){
		if (line[i]==*tab){
			return line.substr(0,i);
		}
	}
	return "";
}
class file_stats{
public:
	vector<int> start_stop;
	map<int, string> relate;
	file_stats(vector<int> Start_stop, map<int, string> Relate){
		start_stop=Start_stop;
		relate=Relate;
	}
};
file_stats getFHstats(string FILE){
	ifstream FH(FILE);
	vector<int> start_stop;
	map<int, string> relate;
	string line;
	string chrom;
	string prev="";
	int i=0; 
	if (FH){
		while (getline(FH,line)){
			chrom = getChrom(line);
			if (prev!=chrom){
				start_stop.push_back(FH.tellg());
				relate[i] 	= chrom;
				i++;
			}
			prev=chrom;
		}

	}else{
		cout<<"\""<<FILE<<"\""<<" doesn't exist, exiting..."<<endl;
	}
	start_stop.push_back(FH.tellg());
	FH.close();
	return file_stats(start_stop, relate);
}
contig::contig(){ this->modified=false; }
void contig::setStats(int st, int sp, double l, double r, double len, float C, string CHROM){
	start 	= st;
	stop 	= sp;
	left 	= l;
	right 	= r;
	cov 	= C ;
	length 	= len;
	chrom 	= CHROM;
    this->modified=true;
}
void contig::display(){
	cout<<chrom<<":"<<start<<"-"<<stop<<endl;
}

vector<double> contig::getVect(bool ChIP){
	vector<double> x;
	if (not ChIP){
        //x.push_back(1);
		//x.push_back(log((left+right)/2));
		//x.push_back(log(length));
		//x.push_back(log(cov/ (left+right+length)));
        
        /*
        if((right+left)/2<=0 || (stop-start+1)<=0 || cov/(stop-start+1)<=0)
        {
            cout<<"Warning: some permutation of values in contig would produce a zero or negative output. l: "<<start<<" r: "<<stop<<" c: "<<cov<<endl;
        }*/
        x.push_back(1);
        x.push_back((start+stop)/2);
        x.push_back((stop-start+1));
        x.push_back(cov/(stop-start+1));
        //x.push_back((right+left)/2);
        //x.push_back((right-left+1));
        //x.push_back(cov/(right-left+1));
	}else{
		x.push_back(1);
		x.push_back(cov / (length+ ((left+right)/2) ));
		x.push_back(cov);
	}
	return x;
}

class contigOut{
public:
	bool EXIT;
	contig * result;
};

contigOut makeContig(string FILE, int start, int stop){
	contigOut CO;
	ifstream FH(FILE);
	FH.seekg(start);
	string line, chrom;
	vector<string> lineArray;
	vector<contig> contigs;
	int st, sp;
	float cov;
	int prevStart= 0;
	int prevStop = 0;
	int p = 0;
	double l = 0;
	double r = 0;
	bool begin=1;
	contig * C;
	C 		= new contig;
	contig * root 	= C;
	float coverage 	= 0;
        int strand=CONTIG_STRAND_POS;

	while (FH.tellg()<stop){
		getline(FH,line);
		lineArray 	= splitter(line, "\t");
		if (lineArray.size()!=4){
                    cout<<endl;
                    cout<<"couldn't parse line: "+line+"\n";
                    cout<<"number of tokens found in given line: "<<lineArray.size()<<endl;
                    CO.EXIT=true;
                    CO.result=NULL;
                    return CO;
		}else{
                    //cout<<line<<endl;

                    chrom 	= lineArray[0];

                    /* NOTE: If you comment out the following block (until the CO.EXIT=true;return CO; nugget)
                    * then FStitch will start to behave like the newer version (ie. it will learn to constantly
                    * switch between states).
                    */
                    //This should be the start position.
                    if(!isNum(lineArray[1]))
                    {
                        cout<<"Element 1 could not be converted"<<endl;
                    }

                    //This should be the stop position
                    if(!isNum(lineArray[2]))
                    {
                        cout<<"Element 2 could not be converted"<<endl;
                    }

                    //This value should represent the frequency of found reads. 
                    if(!isNum(lineArray[3]))
                    {
                        if(stoi(lineArray[3])<0)
                        {
                            //TODO: Consider how to implement better logic for this function.
                            //Perhaps a strand could be specified, then any data not meeting that spec
                            //could then be discarded?
                            strand=CONTIG_STRAND_NEG;
                        }
                        
                        cout<<"Element 3 could not be converted"<<endl;
                    }
                    
                    else
                    {
                        strand=CONTIG_STRAND_POS;
                    }

                    if (not isNum(lineArray[1]) or not  isNum(lineArray[2]) or not isNum(lineArray[3]) ){
                        cout<<endl;
                        cout<<"Line: "<<line<<endl;
                        cout<<"could not convert coordinates or coverage value to number"<<endl;
                        CO.EXIT=true;
                        CO.result=NULL;
                        return CO;
                    }

                    sp 		= stoi(lineArray[2]); //Stop
                    st 		= stoi(lineArray[1]); //Start
                    cov 	= stof(lineArray[3]); //Coverage? Waaait.
                    if (begin){
                        begin 		= 0;
                        prevStop	= sp;
                        l 			= sp-prevStop;
                        prevStart 	= 0;
                    }else if(not begin and (st-p)>2 ){
                        r 			= st - p;
                        C->setStats(prevStart-l,p+r, l, r, p-prevStart, coverage, chrom);
                        C->next 	= new contig;
                        C 			= C->next;
                        prevStart 	= st;
                        l 			= prevStart - p;

                        coverage=0;			
                    }

                    p 		= sp;
                    coverage+=cov;
            }
            
            start++;
	}
	C->next 	= NULL;
	
	FH.close();
	CO.EXIT 	= false;
	CO.result 	= root;
	return CO;
}

contigOut makeContigStrand(string FILE, int start, int stop, int strand){
        contigOut CO;
        ifstream FH(FILE);
        FH.seekg(start);
        string line, chrom;
        vector<string> lineArray;
        vector<contig> contigs;
        int st=0;
        int sp=0;
        float cov=0;
        int prevStart= 0;
        int prevStop = 0;
        int p = 0;
        int i;
        double l = 0;
        double r = 0;
        bool begin=1;
        contig * C=NULL;
        //C->setStats(0, 0, 0, 0, 0, 0, "");
        contig * root;//   = C;
        double coverage  = 0;
        //int strand=CONTIG_STRAND_POS;

        while (FH.tellg()<stop){
                getline(FH,line);
                lineArray       = splitter(line, "\t");
                if (lineArray.size()!=4){
                    cout<<endl;
                    cout<<"couldn't parse line: "+line+"\n";
                    cout<<"Number of tokens found in line: "<<lineArray.size()<<endl;
                    CO.EXIT=true;
                    return CO;
                }else{
                    //cout<<line<<endl;

                    chrom       = lineArray[0];

                    /* NOTE: If you comment out the following block (until the CO.EXIT=true;return CO; nugget)
                    * then FStitch will start to behave like the newer version (ie. it will learn to constantly
                    * switch between states).
                    */
                    //This should be the start position.

                    if(!isNum(lineArray[1]))
                    {
                        cout<<"Element 1 could not be converted"<<endl;
                    }

                    //This should be the stop position
                    if(!isNum(lineArray[2]))
                    {
                        cout<<"Element 2 could not be converted"<<endl;
                    }
                    
                    //cout<<"Value "<<lineArray[3]<<string(isNum(lineArray[3], strand) ? " is" : " is not")<<" a number given strand "<<strand<<endl;
                    
                    if(isNum(lineArray[3], strand))
                    {
                        //cout<<"Added values for line "<<line<<endl;
                        //strand=CONTIG_STRAND_POS;
                        if (not isNum(lineArray[1]) or not  isNum(lineArray[2])){
                            cout<<endl;
                            cout<<"Line: "<<line<<endl;
                            cout<<"could not convert coordinates or coverage value to number given strand "<<strand<<endl;
                            CO.EXIT=true;
                            CO.result=NULL;
                            return CO;
                        }

                        sp          = stoi(lineArray[2]); //Stop
                        st          = stoi(lineArray[1]); //Start
                        cov         = stof(lineArray[3]); //Coverage? Waaait.
                        
                        if(strand)
                        {
                            //If we're negative, then invert the value so fstitch doesn't break.
                            cov*=-1.;
                        }
                        
                        if (begin){
                            begin           = 0;
                            prevStop        = sp;
                            l                       = sp-prevStop;
                            prevStart       = 0;
                        }else if(not begin and (st-p)>2 ){ //NOTE: there are cases in which st is never set.
                            r                       = st - p;
                            //Start, stop, left, right, coverage, length, chromosome.
                            if(!C)
                            {
                                C=new contig();
                                root=C;
                                C->setStats(prevStart-l,p+r, l, r, p-prevStart, coverage, chrom);
                                C->next=NULL;
                            }
                            
                            else
                            {
                                C->next=new contig();
                                C=C->next;
                                C->setStats(prevStart-l,p+r, l, r, p-prevStart, coverage, chrom);
                                C->next=NULL;
                            }
                            
                            
                            //There's a case in which this contig is never properly initialized. 
                            //C->next         = new contig;
                            //C                       = C->next;
                            
                            prevStart       = st;
                            l                       = prevStart - p;

                            coverage=0;                     
                        }

                        p           = sp;
                        coverage+=cov;
                    }
                    
                    //It's clear that not having behavior here is a problem:
                    else
                    {
                        //cout<<"Found line not matching strand."<<endl;
                        //coverage=0;
                    }
            }
            
            start++;
        }
        //Question: why the heck is C null?
        if(C)
        {
            C->next         = NULL;
        }
        
        //This might not necessarily be an undesirable condition: I see in practice that there are loads of contigs with data
        //interspersed with those without.
        else
        {
            cout<<"Potential error condition: C is null. As such, it has never been initialized to a value. Dump of local values and conditions:"<<endl;
            cout<<"Input filename: "<<FILE<<endl;
            cout<<"Start pos: "<<start<<endl;
            cout<<"Stop pos: "<<stop<<endl;
            cout<<"Strand: "<<strand<<endl;
            cout<<"Begin: "<<begin<<endl;
            cout<<"st: "<<st<<endl;
            cout<<"p: "<<p<<endl;
            cout<<"prevStart: "<<prevStart<<endl;
            cout<<"l: "<<l<<endl;
            cout<<"r: "<<r<<endl;
            cout<<"coverage: "<<coverage<<endl;
            cout<<"cov: "<<cov<<endl;
        }
        
        FH.close();
        CO.EXIT         = false;
        CO.result       = root;
        /*
        c=root;
        while(c!=NULL)
        {
            if(!c->modified)
            {
                cout<<"Found unmodified contig. Start: "<<c->start<<" Stop: "<<c->stop<<endl;
                c->start=0;
                c->stop=0;
                c->cov=0;
            }
            
            else
            {
                cout<<"Contig fm "<<c->start<<" to "<<c->stop<<" cov "<<c->cov<<endl; 
            }
            
            c=c->next;
        }*/
        return CO;
}

interval::interval(int st, int sp , string INFO){
	start 	= st;
	stop 	= sp;
	info 	= INFO;
	next 	= NULL;
}
interval::interval(){}

readTrainingFileReturn::readTrainingFileReturn(){}

readTrainingFileReturn readTrainingFile(string FILE){
	readTrainingFileReturn RETURN;
	map<string, interval *> 	R;
	map<string, interval *> 	roots;
	
	vector<int> start_stop(3);
	ifstream FH(FILE);
	string line;
	vector<string>lineArray;
	if (FH){
		while (getline(FH, line)){
			lineArray 			= splitter(line, "\t");
			if (lineArray.size()!=4){
				cout<<"Line: "<<line<<", in training file is not formatted properly, tab delimited"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}

			if (not (lineArray[3]=="0" or lineArray[3]=="1")){
				cout<<endl;
				cout<<"Line: "<<line<<", must have training label as either 0 or 1"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}
			if (not isNum(lineArray[1]) or not isNum(lineArray[2])){
				cout<<"Line: "<<line<<", coordinates must be numbers"<<endl;
				RETURN.EXIT 	= true;
				return RETURN;
			}
			if (R.find(lineArray[0])==R.end()){
				R[lineArray[0]] 		= new interval(stoi(lineArray[1] ), stoi(lineArray[2]), lineArray[3]);
				roots[lineArray[0]] 	= R[lineArray[0]];
			}else{
				R[lineArray[0]]->next 	= new interval(stoi(lineArray[1] ), stoi(lineArray[2]), lineArray[3]);
				R[lineArray[0]] 		= R[lineArray[0]]->next;
			}

		}
	}else{
		cout<<"couldn't open: "<<FILE<<endl;
		cout<<"exiting..."<<endl;
	}
	typedef map<string, interval *>::iterator r_it;
	for (r_it I = R.begin(); I!=R.end(); I++){
		R[I->first] 	= roots[I->first];
	}

	FH.close();
	RETURN.EXIT 	= false;
	RETURN.result 	= R;
	return RETURN;

}

int lineCompare(string line1, string line2)
{
    char *line1c, *line2c;
    char *ptr;
    int l1col, l2col;
    
    line1c=(char*) line1.c_str();
    line2c=(char*) line2.c_str();

    //Awful kluge to find the second column in each line:
    for(ptr=line1c;(*ptr)!='\t';ptr++);
    ptr++;
    l1col=atoi(ptr);

    for(ptr=line2c;(*ptr)!='\t';ptr++);
    ptr++;
    l2col=atoi(ptr);

    return l1col-l2col;
}

int histogramLineCompare(string line1, string line2, string prevWrittenLine)
{
    int cmpV;
    vector<string> l1toks=splitter(line1, "\t");
    vector<string> l2toks=splitter(line2, "\t");
    
    //First, determine if these are for the same chromosome:
    if(l1toks[0]==l2toks[0])
    {
        return atoi(l1toks[1].c_str())-atoi(l2toks[1].c_str());
    }
    
    //If they aren't, then give priority to the most similar line.
    else
    {
        vector<string> prevLineToks=splitter(prevWrittenLine, "\t");
        
        //A negative value indicates that line 1 has a lesser starting position than line 2 and should be written first.
        if(prevLineToks[0]==l1toks[0])
        {
            return -1;
        }
        
        //A positive value indicates that line 1 has a greater starting position than line 2 and should be written last.
        else
        {
            return 1;
        }
    }
    
}

bool fileExists(string filename)
{
    struct stat s;
    
    return stat(filename.c_str(), &s)==0;
}

readTrainingFileReturn readSplitTrainingFile(string onfile, string offfile){
        readTrainingFileReturn RETURN;
        map<string, interval *>         R;
        map<string, interval *>         roots;

        vector<int> start_stop(3);
        //ifstream FH(FILE);
        ifstream onFH(onfile);
        ifstream offFH(offfile);
        string line;
        vector<string> lineArray;
        vector<string> lines;
        bool ongood;
        bool offgood;
        string online;
        string offline;
        int cmpV;
        int i;

        // Fill up the vector of lines while taking advantage of the fact that the input
        // files *should* be sorted.
        if(onFH && offFH)
        {
            ongood=getline(onFH, online).good();
            offgood=getline(offFH, offline).good();
            
            while(ongood&&offgood)
            {
                cmpV=lineCompare(online, offline);
                if(cmpV>0) //Ie. the starting position of line1>line2.
                {
                    lines.push_back(offline+"\t0");
                    offgood=getline(offFH, offline).good();
                }

                //Ie. the starting position of line1=line2.
                else if(cmpV==0)
                {
                    lines.push_back(offline+"\t0");
                    lines.push_back(online+"\t1");
                    offgood=getline(offFH, offline).good();
                    ongood=getline(onFH, online).good();
                }

                //Ie. the starting position of line1<line2.
                else
                {
                    lines.push_back(online+"\t1");
                    ongood=getline(onFH, online).good();
                }
            }

            //Attempt to push the rest of the data into the vector:
            while(ongood)
            {
                lines.push_back(online+"\t1");
                ongood=getline(onFH, online).good();
            }

            while(offgood)
            {
                lines.push_back(offline+"\t0");
                offgood=getline(offFH, offline).good();
            }
            
            cout << "Got " << lines.size() << " lines in split training file." << endl;

            //Now that we have a proper set of lines:
            for(i=0;i<lines.size();i++)
            {
                line=lines[i];

                //Code shamelessly copy-pasted from above:

                if (!line.empty() && line[line.size() - 1] == '\r'){
                        line.erase(line.size() - 1);
                }
                lineArray                       = splitter(line, "\t");
                if (lineArray.size()!=4){
                        cout<<"Line: "<<line<<", in training file is not formatted properly, tab delimited"<<endl;
                        RETURN.EXIT     = true;
                        return RETURN;
                }
                //NOTE: The splitter class doesn't remove carriage returns before line breaks:
                if (not (lineArray[3].at(0)=='0' or lineArray[3].at(0)=='1')){
                        printf("Line: %s, must contain either 0 or 1 as training input\n",line.c_str() );
                        RETURN.EXIT     = true;
                        return RETURN;
                }
                if (not isNum(lineArray[1]) or not isNum(lineArray[2])){
                        cout<<"Line: "<<line<<", coordinates must be numbers"<<endl;
                        RETURN.EXIT     = true;
                        return RETURN;
                }
                if (R.find(lineArray[0])==R.end()){
                        R[lineArray[0]]                 = new interval(stoi(lineArray[1] ), stoi(lineArray[2]), lineArray[3]);
                        roots[lineArray[0]]     = R[lineArray[0]];
                }else{
                        R[lineArray[0]]->next   = new interval(stoi(lineArray[1] ), stoi(lineArray[2]), lineArray[3]);
                        R[lineArray[0]]                 = R[lineArray[0]]->next;
                }
            }
        }

        else
        {
            cerr<<"Couldn't open either "<<onfile<<" or "<<offfile<<"."<<endl;
        }

        typedef map<string, interval *>::iterator r_it;
        for (r_it I = R.begin(); I!=R.end(); I++){
                R[I->first]     = roots[I->first];
        }

        onFH.close();
        offFH.close();
        RETURN.EXIT     = false;
        RETURN.result   = R;
        return RETURN;
}

map<string,contig *> readBedGraphFile(string FILE, map<string, interval *> T, bool verbose){
	file_stats fs 					= getFHstats(FILE);
	vector<int> start_stop 			= fs.start_stop;

	map<int,string> relate 			= fs.relate;
	map<string,contig *> 	D;
	vector<contig *> M(start_stop.size());
	int nthreads 					= omp_get_max_threads();
	bool abort = false;
	#pragma omp parallel for
	for(int n=1; n<start_stop.size(); ++n)
	{
		#pragma omp flush (abort)
		if (T.find(relate[n-1])!=T.end()){
			contigOut CO 	= makeContig(FILE, start_stop[n-1],start_stop[n]);
			if (CO.EXIT){
				D.clear();
				abort = true;
			}
			M[n-1] 			= CO.result;
			D[relate[n-1]] 	= M[n-1];
		}
		
	}
	if (abort){
		D.clear();
	}
	return D;
}

map<string,contig *> readBedGraphFileStrand(string FILE, map<string, interval *> T, bool verbose, int strand){
        file_stats fs                                   = getFHstats(FILE);
        vector<int> start_stop                  = fs.start_stop;

        map<int,string> relate                  = fs.relate;
        map<string,contig *>    D;
        vector<contig *> M(start_stop.size());
        int nthreads                                    = omp_get_max_threads();
        bool abort = false;
        #pragma omp parallel for
        for(int n=1; n<start_stop.size(); ++n)
        {
                #pragma omp flush (abort)
                if (T.find(relate[n-1])!=T.end()){
                        contigOut CO    = makeContigStrand(FILE, start_stop[n-1],start_stop[n], strand);
                        if (CO.EXIT){
                                D.clear();
                                abort = true;
                        }
                        M[n-1]                  = CO.result;
                        D[relate[n-1]]  = M[n-1];
                }
                
        }
        if (abort){
                D.clear();
        }
        return D;
}

//NOTE: We should be sorting by chromosome first, then position within the chromosome.
tmpfile_t mergeSplitBedGraph(string posFile, string negFile)
{
    tmpfile_t retTmpFile;
    int fd;
    ifstream posInFile(posFile);
    ifstream negInFile(negFile);
    bool posGood;
    bool negGood;
    string posLine;
    string negLine;
    int cmpV;
    string prevLine;
    
    retTmpFile.tempName=strdup("tmpXXXXXX");
    fd=mkstemp(retTmpFile.tempName);
    
    retTmpFile.tempFile=fdopen(fd, "w+");
    
    if(posInFile && negInFile)
    {
        posGood=getline(posInFile, posLine).good();
        negGood=getline(negInFile, negLine).good();
        
        //For the first line:
        prevLine=posGood;
        
        while(posGood && negGood)
        {
            cmpV=histogramLineCompare(posLine, negLine, prevLine);
            if(cmpV>0) //Ie. the starting position of line1>line2.
            {
                fprintf(retTmpFile.tempFile, "%s\n", negLine.c_str());
                prevLine=negLine;
                negGood=getline(negInFile, negLine).good();
            }

            //Ie. the starting position of line1=line2.
            else if(cmpV==0)
            {
                fprintf(retTmpFile.tempFile, "%s\n", negLine.c_str());
                fprintf(retTmpFile.tempFile, "%s\n", posLine.c_str());
                prevLine=posLine;
                posGood=getline(posInFile, posLine).good();
                negGood=getline(negInFile, negLine).good();
            }

            //Ie. the starting position of line1<line2.
            else
            {
                fprintf(retTmpFile.tempFile, "%s\n", posLine.c_str());
                prevLine=posLine;
                posGood=getline(posInFile, posLine).good();
            }
        }

        //Attempt to push the rest of the data into the vector:
        while(posGood)
        {
            fprintf(retTmpFile.tempFile, "%s\n", posLine.c_str());
            posGood=getline(posInFile, posLine).good();
        }

        while(negGood)
        {
            fprintf(retTmpFile.tempFile, "%s\n", negLine.c_str());
            negGood=getline(negInFile, negLine).good();
        }
        
        fflush(retTmpFile.tempFile);
    }
    
    else
    {
        fclose(retTmpFile.tempFile);
        retTmpFile.tempFile=NULL;
    }
        
    return retTmpFile;
}

vector<string> readFileLines(string fileName)
{
    ifstream inF(fileName);
    string str;
    vector<string> outVect;
    
    if(inF)
    {
        while(getline(inF, str).good())
        {
            outVect.push_back(str);
        }
    }
    
    return outVect;
}

map<string,contig *> readSplitBedGraphFileStrand(string posFile, string negFile, map<string, interval *> T, bool verbose, int strand)
{
    tmpfile_t t;
    
    t=mergeSplitBedGraph(posFile, negFile);
    
    //The merge was unsuccessful. 
    if(!t.tempFile)
    {
        free(t.tempName);
        map<string, contig*> r;
        cout<<"Error: merge unsuccessful."<<endl;
        return r;
    }
    
    else
    {
        string FILE(t.tempName);
        map<string, contig*> r=readBedGraphFileStrand(FILE, T, verbose, strand);
        fclose(t.tempFile);
        remove(t.tempName);
        cout<<"Removing temporary merge file: "<<t.tempName<<endl;
        free(t.tempName);
        
        return r;
    }
}

map<string,contig *> readBedGraphFileAll(string FILE,int np){
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(np); // Use 4 threads for all consecutive parallel regions
	
	file_stats fs 					= getFHstats(FILE);
	vector<int> start_stop 			= fs.start_stop;
	map<int,string> relate 			= fs.relate;
	map<string,contig *> 	D;
	vector<contig *> M(start_stop.size());
	bool abort = false;
	#pragma omp parallel for
	for(int n=1; n<start_stop.size(); ++n)
	{
		#pragma omp flush (abort)
		contigOut CO 	= makeContig(FILE, start_stop[n-1],start_stop[n]);
		if (CO.EXIT){
			D.clear();
			abort = true;
		}
		M[n-1] 			= CO.result;
		D[relate[n-1]] 	= M[n-1];
		
	}
	if (abort){
		D.clear();
	}
	return D;
}

//This follows the conventions for strand established in main_segment:
//"+" for pos strand, "-" for neg strand, and "." for both (unsupported).
map<string, contig *> readBedGraphFileAllGivenStrand(string FILE, int np, string strand)
{
    int str;
    
    str=(strand!="+");
    
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(np); // Use 4 threads for all consecutive parallel regions
	
	file_stats fs 					= getFHstats(FILE);
	vector<int> start_stop 			= fs.start_stop;
	map<int,string> relate 			= fs.relate;
	map<string,contig *> 	D;
	vector<contig *> M(start_stop.size());
	bool abort = false;
	#pragma omp parallel for
	for(int n=1; n<start_stop.size(); ++n)
	{
		#pragma omp flush (abort)
		contigOut CO 	= makeContigStrand(FILE, start_stop[n-1], start_stop[n], str);//makeContig(FILE, start_stop[n-1],start_stop[n]);
		if (CO.EXIT){
			D.clear();
			abort = true;
		}
		M[n-1] 			= CO.result;
		D[relate[n-1]] 	= M[n-1];
		
	}
	if (abort){
		D.clear();
	}
	return D;
}


map<string, map<string, interval *>> readRefSeq(string FILE){
	map<string, map<string, interval *>> R;	
	map<string, map<string, interval *>> roots;	
	
	R["+"] =	map<string, interval *>();
	R["-"] = 	map<string, interval *>();
	roots["+"] =	map<string, interval *>();
	roots["-"] = 	map<string, interval *>();
	
	if (not FILE.empty()){
		ifstream FH(FILE);
		string line, chrom, strand, ID;
		int 	start, stop; 
		vector<string> lineArray;
		if (FH){
			while (getline(FH,line)){
				lineArray 	= splitter(line, "\t");
				chrom 		= lineArray[2], strand = lineArray[3], ID=lineArray[1];
				start 		= stoi(lineArray[4]), stop=stoi(lineArray[5]);
				if (R[strand].find(chrom)==R[strand].end()){
					R[strand][chrom] 		= new interval(start, stop, ID);
					roots[strand][chrom] 	= R[strand][chrom];
				}else{
					R[strand][chrom]->next 	= new interval(start, stop, ID);
					R[strand][chrom] 		= R[strand][chrom]->next;
				}
			}
		}else{
			cout<<"couldn't open: "<<FILE<<"\nExiting..."<<endl;
			return R;
		}
		FH.close();
	}
	typedef map<string, map<string, interval *>>::iterator st_it;
	typedef map<string, interval *>::iterator chrom_it;
	
	for (st_it strand = R.begin(); strand != R.end(); strand++ ){
		for (chrom_it chrom = strand->second.begin(); chrom!= strand->second.end(); chrom++){
			R[strand->first][chrom->first] 	= roots[strand->first][chrom->first];
		}
	}
	
	return R;
}
RTOF::RTOF(vector<double> w, vector<vector<double>> a, bool CH){
		W=w,A=a;
		EXIT=false;
		ChIP=CH;
}
RTOF::RTOF(){
	EXIT=true;
}

RTOF readTrainingOutFile(string FILE){
	ifstream FH(FILE);
	string line;
	vector<double> W;
	vector<vector<double>> A;
	vector<string> lineArray;
        string commandLine="not defined";
	bool begin 	= 1;
	bool ChIP 	= 0;
        RTOF retRTOF;
	if (FH){
		while (getline(FH,line)){
			if (begin and ("#" != line.substr(0,1) ) ){
				RTOF ROOT;
				cout<<"This is not an output training file\nfrom the fast read stitcher"<<endl;
				return ROOT;
			}

			begin = false;
			if ("#ChIP" == line.substr(0,5)){
				lineArray 		= splitter(line, ":");
				ChIP 			= (lineArray[1]=="1");
			}
			
                        // Read the training command line and put it in the output struct.
			else if(line.substr(0,8)=="#Command")
                        {
                            lineArray=splitter(line, ":");
                            commandLine=lineArray[1];
                        }
                        
			if ("#" != line.substr(0,1)){
				lineArray 		= splitter(line, ":");
				if (lineArray[0].substr(0,1)=="L"){
					lineArray 	= splitter(lineArray[1], ",");
					for (int i = 0; i < lineArray.size(); i++){
						W.push_back(stof(lineArray[i]));
					}
				}else if(lineArray[0].substr(0,1)=="H"){
					lineArray 	= splitter(lineArray[1], ",");
					int k 		= 0;
					for (int i  = 0; i < 2;i++){
						vector<double> row;
						for (int j=0;j<2;j++){
							row.push_back(stof(lineArray[k]));
							k++;
						}
						A.push_back(row);
					}
				}
			}
		}

	}else{
		cout<<"\""<<FILE<<"\""<<" doesn't exist, exiting..."<<endl;
	}
	FH.close();
        
        retRTOF=RTOF(W, A, ChIP);
        retRTOF.commandLine=commandLine;
	return retRTOF;
}

int getStrand(string line)
{
    double val;
    vector<string> lineToks;
    
    lineToks=splitter(line, "\t");
    
    if(lineToks.size()<4)
    {
        return 0;
    }
    
    else
    {
        val=strtod(lineToks[3].c_str(), NULL);
        
        if(val<0)
        {
            return -1;
        }
        
        else if(val>0)
        {
            return 1;
        }
        
        else
        {
            return 0;
        }
    }
}

/* This function is a last resort used by main_segment and main_train to determine whether or not a given input histogram contains positive, negative, or mixed positive and negative reads. 
 * 
 * Return values:
 * STRAND_UNSPECIFIED on error
 * STRAND_POSITIVE on positive
 * STRAND_NEGATIVE on negative
 * STRAND_BOTH on both.
 */
int checkBedFileType(string bedName)
{
    int firstLineStrand;
    ifstream f(bedName);
    int strand;
    string l;
    
    if(f)
    {
        getline(f, l);
        firstLineStrand=getStrand(l);
        
        //If our first line has 0 coverage, then we need to keep fetching lines until we find something with coverage.
        while(!firstLineStrand && getline(f, l))
        {
            firstLineStrand=getStrand(l);
        }
        
        while(getline(f, l))
        {
            strand=getStrand(l);
            
            if(strand!=firstLineStrand && strand)
            {
                return STRAND_BOTH;
            }
        }
        
        if(firstLineStrand<0)
        {
            return STRAND_NEGATIVE;
        }
        
        else if(firstLineStrand>0)
        {
            return STRAND_POSITIVE;
        }
        
        else
        {
            return STRAND_UNSPECIFIED;
        }
    }
    
    else
    {
        return STRAND_UNSPECIFIED;
    }
}




