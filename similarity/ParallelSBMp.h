#ifndef PARALLELSBMP_H
#define PARALLELSBMP_H

#include <vector>
#include <omp.h>
#include "AscGrid.h"
#include "AttributeRule.h"
#include "CategoryIntegration.h"
#include "SampleIntegration.h"

#define VERY_SMALL 0.000001
using namespace std;

class ParallelSBMp
{
private:

	vector<AttributeRule> attriRules;
	CategoryIntegration catIntegration;
	SampleIntegration sampleIntegration;
	double uncertaintyThreshold;

	double * sampleClimateV;
	double * sampleGeologyV;
	double * sampleTerrainV;
	double * sampleVegetationV;
	double * sampleOtherV;

	double * climateVRange;
	double * geologyVRange;
	double * terrainVRange;
	double * vegeVRange;
	double * otherVRange;

	int rowIndex;
	int colIndex;
	double curX;
	double curY;
	double lowerLeftX;
	double lowerLeftY;
	double cellSize;
	int totalRows;

	double string_to_double( const std::string& s );

public:
	ParallelSBMp();
	ParallelSBMp(double * sampleClimateV, double * sampleGeologyV, double * sampleTerrainV,
		double * sampleVegetationV, double * sampleOtherV, vector<string> &AttriRules,
		double * climateVRange,double * geologyVRange,double * terrainVRange,double * vegeVRange,double * otherVRange,
		int rowIndex, int colIndex, double curX, double curY, double lowerLeftX, double lowerLeftY, double cellSize, int totalRows,
		string catIntegrationMethod, string sampleIntegrationMethod, string uncertaintyThreshold);

	void getPropertyMap(double *climateStd,double *geoStd,double *terrainStd,double *vegetationStd,double *otherStd,
		double ** climate2Dvalues, double ** geo2Dvalues, double ** terrain2Dvalues, double ** vege2Dvalues, double ** other2Dvalues,
		double ** similarityVals, double ** uncertaintyVals,
		int ClimateLyrCnt, int GeogLyrCnt, int TerrainLyrCnt, int VegeLyrCnt, int OtherLyrCnt,int rows, int cols,
		double noData, double *climateNoData,double *geoNoData,double *terrainNoData,double *vegetationNoData,double *otherNoData,
		double *climateSum, double *geologySum, double *terrainSum, double *vegetationSum, double *otherSum,
		int *climateCount, int *geologyCount, int *terrainCount, int *vegetationCount, int *otherCount);//wtf
	void split(string s, string sep, vector<string> &flds);
};

#endif
