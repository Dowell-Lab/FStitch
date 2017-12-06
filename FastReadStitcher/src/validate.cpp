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
    int val;
    val=strtol(num.c_str(), NULL, 10);
    
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
