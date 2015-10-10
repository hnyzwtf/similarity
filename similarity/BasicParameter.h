#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include "AscGrid.h"
#include "AttributeRule.h"
#include "ParallelSBMp.h"
#include <sstream>
#include <fstream>
#include <iomanip>
using namespace std;
class BasicParameter
{
public:
	BasicParameter();
	BasicParameter(char *attriRules, char *environLyrsPath, char *sample, char *uncertaintyThreshold, char *simiLayerPath, char *uncertaintyLayerPath, 
		const char *categoryIntegrationMethod, const char *sampleIntegrationMethod);
	void getBasicParameter();
	void parseStr(string str, char c, vector<string>& tokens);
	double string_to_double( const std::string& s );
	static double string_to_double1( const std::string& s );
	double *getSum(double **category2DValues, int categoryLyrCnt, double *categoryNoData, double *sampleCategoryV, int rows, int cols);
	int *getCount(double **category2DValues, int categoryLyrCnt, double *categoryNoData, int rows, int cols);

private:
	char *attriRules;
	char *environLyrsPath;
	char *sample;
	char *uncertaintyThreshold;
	char *simiLayerPath;
	char *uncertaintyLayerPath;
	const char *categoryIntegrationMethod;
	const char *sampleIntegrationMethod;
	
};