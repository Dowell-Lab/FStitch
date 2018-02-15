#ifndef read_H
#define read_H
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <sys/stat.h>

/* For strand definitions. */
#include "ParamWrapper.h"
using namespace std;

#define CONTIG_STRAND_POS 1
#define CONTIG_STRAND_NEG 2
class contig{
public:
	int start, stop;
	double left, right, length;
	float cov;
	string chrom;
	contig * next;
	contig();
	void setStats(int , int , double , double , double, float , string);
	void display();
	vector<double> getVect(bool);
        int strand;
    bool modified;
};
class interval{
public:
	int start;
	int stop;
	string info;
	interval * next = NULL;
	interval(int, int, string);
	interval();	
};
class RTOF{
public:
	vector<double> W;
	vector<vector<double>> A;
	bool ChIP;
	bool EXIT;
	RTOF(vector<double>,vector<vector<double>>, bool);
	RTOF();
        string commandLine;
};

class readTrainingFileReturn{
public:
	bool EXIT;
	map<string, interval *> result;
	readTrainingFileReturn();
	
};

typedef struct
{
    string in1;
    char in1ch[80];
    string in2;
    char in2ch[80];
    
    bool splitGenerated=false;

    bool wasSplit=false;
} splitinputfile_t;

typedef struct
{
    char *tempName;
    FILE *tempFile;
} tmpfile_t;

splitinputfile_t splitBedgraphFile(string);
bool fileExists(string filename);
readTrainingFileReturn readTrainingFile(string);
int lineCompare(string line1, string line2);
int histogramLineCompare(string line1, string line2, string prevWrittenLine);
readTrainingFileReturn readSplitTrainingFile(string onfile, string offfile);
map<string,contig *> readBedGraphFile(string, map<string, interval *>, bool);
map<string,contig *> readBedGraphFileStrand(string, map<string, interval *>, bool, int);
tmpfile_t mergeSplitBedGraph(string, string);
map<string,contig *> readSplitBedGraphFileStrand(string, string, map<string, interval *>, bool, int);
map<string, map<string, interval *>> readRefSeq(string);
RTOF readTrainingOutFile(string);
map<string, contig *> readBedGraphFileAllGivenStrand(string FILE, int np, string strand);
map<string,contig *> readBedGraphFileAll(string,int);
vector<string> readFileLines(string fileName);

/* This function determines which strand is represented by a given histogram bedgraph line.
 * In other words, it determines if the last field in line, when tokenized, is positive or negative.
 * 
 * This function returns 1 on positive, -1 on negative.
 */
int getStrand(string line);

/* This function determines whether a given histogram bedgraph file contains positive reads, 
 * negative reads, or both positive and negative reads.
 * 
 * This function returns STRAND_POSITIVE for a positive strand, STRAND_NEGATIVE for a negative strand,
 * STRAND_BOTH for both positive and negative strands, and STRAND_UNSPECIFIED if there is an error in reading
 * or parsing the given file.
 */
int checkBedFileType(string bedName);

#endif
