#ifndef read_H
#define read_H
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <sys/stat.h>
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

    bool wasSplit=false;
} splitinputfile_t;

splitinputfile_t splitBedgraphFile(string);
bool fileExists(string filename);
readTrainingFileReturn readTrainingFile(string);
int lineCompare(string line1, string line2);
readTrainingFileReturn readSplitTrainingFile(string onfile, string offfile);
map<string,contig *> readBedGraphFile(string, map<string, interval *>, bool);
map<string,contig *> readBedGraphFileStrand(string, map<string, interval *>, bool, int);
map<string, map<string, interval *>> readRefSeq(string);
RTOF readTrainingOutFile(string);
map<string,contig *> readBedGraphFileAll(string,int);

#endif
