#include "validate.h"
#include <string>
#include <locale>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

bool isNum(string num){
	for (int i=0; i < num.size(); i++){
		if ((not isdigit(num[i])) and num.substr(i,1) != "." ){
			return false;
		}
	}
	return true;
}

bool isNum(string num, int strand)
{ 
    //So we're totally not ambiguous at all:
    long val;
    //strtol doesn't work with scientific notation. This is just great...
    val=(long) strtod(num.c_str(), NULL);//strtol(num.c_str(), NULL, 10);
    
    if(errno==EINVAL||errno==ERANGE)
    {
        return false;
    }
    
    //If we're looking at the negative strand:
    if(strand)
    {
        return val<=0;
    }
    
    else
    {
        return val>=0;
    }
}

bool isFile(string FILE){
	ifstream f(FILE);
	if (f.good()) {
		f.close();
		return true;
	}else{
		f.close();
		return false;
	}   
}
