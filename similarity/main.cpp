
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include "BasicParameter.h"
#include <sstream>
#include <fstream>
#include <iomanip>
using namespace std;

int main(int argc, char **argv)
{
	GDALAllRegister();	
	 //the input parameter
	char *environLyrsPath = argv[1];
	char *sample = argv[2]; 
	char * attriRules = argv[3];
	char * threshold = argv[4];
	char * simiPath = argv[5];
	char * unLayerPath = argv[6];
		
    char *uncertaintyThreshold = threshold;
    char *simiLayerPath = simiPath;
    char *uncertaintyLayerPath = unLayerPath;
	const char *categoryIntegrationMethod = "Limit";
	const char *sampleIntegrationMethod = "Limit";
	
	BasicParameter inf(attriRules, environLyrsPath, sample, uncertaintyThreshold, simiLayerPath,
		uncertaintyLayerPath, categoryIntegrationMethod,sampleIntegrationMethod);
	inf.getBasicParameter();//wtf
	//getchar();
	return 0;

}
