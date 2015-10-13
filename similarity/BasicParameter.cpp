#include "BasicParameter.h"
#include <sstream>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
using namespace std;

void BasicParameter::parseStr(string str, char c, vector<string>& tokens){
	unsigned int posL = 0;
	unsigned int posR = 0;
	while(posR < str.length()-1){
		posR = str.find_first_of(c,posL);
		string sub = str.substr(posL,posR-posL);
		tokens.push_back(sub);
		posL = posR + 1;
	}
}
double BasicParameter:: string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

double * BasicParameter::getSum(double **category2DValues, int categoryLyrCnt, double *categoryNoData, double *sampleCategoryV, int rows, int cols)
{
	double *sum = NULL;
	sum = new double[categoryLyrCnt];
	for (int i = 0; i <categoryLyrCnt; i++)
	{
		sum[i] = 0;
		for (int rowId = 0; rowId < rows; rowId++)
		{
			for (int colId = 0; colId < cols; colId++)
			{
				double gridV = category2DValues[i][rowId * cols + colId];
				if (abs(gridV - categoryNoData[i]) >= VERY_SMALL)
				{
					sum[i] += pow(gridV - sampleCategoryV[i], 2);
				}
			}
		}
	}
	return sum;
}
int * BasicParameter::getCount(double **category2DValues, int categoryLyrCnt, double *categoryNoData, int rows, int cols)
{
	int *categoryCounter = NULL;
	categoryCounter = new int[categoryLyrCnt];
	for (int i = 0; i <categoryLyrCnt; i++)
	{
		categoryCounter[i] = 0;
		for (int rowId = 0; rowId < rows; rowId++)
		{
			for (int colId = 0; colId < cols; colId++)
			{
				double gridV = category2DValues[i][rowId * cols + colId];
				if (abs(gridV - categoryNoData[i]) >= VERY_SMALL)
				{
					categoryCounter[i] += 1;
				}
			}
		}
	}
	return categoryCounter;
}

BasicParameter::BasicParameter()
{
	
}
BasicParameter::BasicParameter(char *attriRules, char *environLyrsPath, char *sample, char *uncertaintyThreshold, char *simiLayerPath, 
	char *uncertaintyLayerPath, const char *categoryIntegrationMethod, const char *sampleIntegrationMethod)
{
	this->attriRules = attriRules;
	this->environLyrsPath = environLyrsPath;
	this->sample = sample;
	this->uncertaintyThreshold = uncertaintyThreshold;
	this->simiLayerPath = simiLayerPath;
	this->uncertaintyLayerPath = uncertaintyLayerPath;
	this->categoryIntegrationMethod = categoryIntegrationMethod;
	this->sampleIntegrationMethod = sampleIntegrationMethod;
	
}

void BasicParameter::getBasicParameter()
{
	
	// the count of input layer
	int climateLyrCnt = 0;
	int geologyLyrCnt = 0;
	int terrainLyrCnt = 0;
	int vegeLyrCnt = 0;
	int otherLyrCnt = 0;
	int totalRows = 0;
	int totalCols = 0;
	// attribute rules
	vector<string> attributeRules;// 以#分割我们输入的规则
	string attriRulesList = string(attriRules);
	parseStr(attriRulesList,'#',attributeRules);
	vector<string> typeOfLayers;
	for(int i=0; i<attributeRules.size(); i++)
	{
		vector<string> ruleParameterArray;
		parseStr(attributeRules[i],'?',ruleParameterArray);
		string type = ruleParameterArray[0];
		typeOfLayers.push_back(type);
		if(type == "Climate")
			climateLyrCnt++;
		else if(type == "Geology")
			geologyLyrCnt++;
		else if(type == "Terrain")
			terrainLyrCnt++;
		else if(type == "Vegetation")
			vegeLyrCnt++;
		else if(type == "Other")
			otherLyrCnt++;
	}
	// environment layers
	AscGrid *climateLyrs = NULL;
	AscGrid *geologyLyrs = NULL;
	AscGrid *terrainLyrs = NULL;
	AscGrid *vegetationLyrs = NULL;
	AscGrid *otherLyrs = NULL;
	// range
	double *climateVRange = NULL;
	double *geologyVRange = NULL;
	double *terrainVRange = NULL;
	double *vegeVRange = NULL;
	double *otherVRange = NULL;
	// the cellvalue of the layer in the position of sample index
	double *sampleClimateV = NULL;
	double *sampleGeologyV = NULL;
	double *sampleTerrainV = NULL;
	double *sampleVegeV = NULL;
	double *sampleOtherV = NULL;
	
	double noData = -9999;
	double lowerLeftY;
	double lowerLeftX;
	double cellSize;

	int climateCnt = 0;
	int geologyCnt = 0;
	int terrainCnt = 0;
	int vegeCnt = 0;
	int otherCnt = 0;

	double *climateNoData = NULL;
	double *geologyNoData = NULL;
	double *terrainNoData = NULL;
	double *vegeNoData = NULL;
	double *otherNoData = NULL;

	double *climateStd = NULL;
	double *geologyStd = NULL;
	double *terrainStd = NULL;
	double *vegeStd = NULL;
	double *otherStd = NULL;

	double **climate2DValues = NULL;
	double **geology2DValues = NULL;
	double **terrain2DValues = NULL;
	double **vege2DValues = NULL;
	double **other2DValues = NULL;

	// load all environment layers
	vector<string> environLyrs;
	parseStr(string(environLyrsPath), '#', environLyrs);
	int numOfLayers = environLyrs.size();
	if (climateLyrCnt > 0)
	{
		climateLyrs = new AscGrid[climateLyrCnt];
		climate2DValues = new double *[climateLyrCnt];
		climateVRange = new double[climateLyrCnt];
		climateNoData = new double[climateLyrCnt];
		climateStd = new double[climateLyrCnt];
	}
	if (geologyLyrCnt > 0)
	{
		geologyLyrs = new AscGrid[geologyLyrCnt];
		geology2DValues = new double *[geologyLyrCnt];
		geologyVRange = new double[geologyLyrCnt];
		geologyNoData = new double[geologyLyrCnt];
		geologyStd = new double[geologyLyrCnt];
	}
	if (terrainLyrCnt > 0)
	{
		terrainLyrs = new AscGrid[terrainLyrCnt];
		terrain2DValues = new double *[terrainLyrCnt];
		terrainVRange = new double[terrainLyrCnt];
		terrainNoData = new double[terrainLyrCnt];
		terrainStd = new double[terrainLyrCnt];
	}
	if (vegeLyrCnt > 0)
	{
		vegetationLyrs = new AscGrid[vegeLyrCnt];
		vege2DValues = new double *[vegeLyrCnt];
		vegeVRange = new double[vegeLyrCnt];
		vegeNoData = new double[vegeLyrCnt];
		vegeStd = new double[vegeLyrCnt];
	}
	if (otherLyrCnt > 0)
	{
		otherLyrs = new AscGrid[otherLyrCnt];
		other2DValues = new double *[otherLyrCnt];
		otherVRange = new double[otherLyrCnt];
		otherNoData = new double[otherLyrCnt];
		otherStd = new double[otherLyrCnt];
	}
	for (int i = 0; i < numOfLayers; i++)
	{
		AscGrid attributeLayer;
		attributeLayer.readAscGridGDAL(environLyrs[i]);
		string type = typeOfLayers[i];
		double range = attributeLayer.getMax() - attributeLayer.getMin();
		if (type == "Climate")
		{
			climateLyrs[climateCnt] = attributeLayer;
			climateVRange[climateCnt] = range;
			climateNoData[climateCnt] = attributeLayer.getNodaVal();
			climateStd[climateCnt] = attributeLayer.getStdValue();
			climateCnt++;
		}
		else if (type == "Geology")
		{
			geologyLyrs[geologyCnt] = attributeLayer;
			geologyVRange[geologyCnt] = range;
			geologyNoData[geologyCnt] = attributeLayer.getNodaVal();
			geologyStd[geologyCnt] = attributeLayer.getStdValue();
			geologyCnt++;
		}
		else if (type == "Terrain")
		{
			terrainLyrs[terrainCnt] = attributeLayer;
			terrainVRange[terrainCnt] = range;
			terrainNoData[terrainCnt] = attributeLayer.getNodaVal();
			terrainStd[terrainCnt] = attributeLayer.getStdValue();
			terrainCnt++;
		}
		else if (type == "Vegetation")
		{
			vegetationLyrs[vegeCnt] = attributeLayer;
			vegeVRange[vegeCnt] = range;
			vegeNoData[vegeCnt] = attributeLayer.getNodaVal();
			vegeStd[vegeCnt] = attributeLayer.getStdValue();
			vegeCnt++;
		}
		else if (type == "Other")
		{
			otherLyrs[otherCnt] = attributeLayer;
			otherVRange[otherCnt] = range;
			otherNoData[otherCnt] = attributeLayer.getNodaVal();
			otherStd[otherCnt] = attributeLayer.getStdValue();
			otherCnt++;
		}
	}  // end of for (int i = 0; i < numOfLayers; i++)
	// compute some parameters in each raster
	if (climateLyrCnt > 0)
	{
		totalRows = climateLyrs[0].getNumOfRows();
		totalCols = climateLyrs[0].getNumOfCols();
		lowerLeftX = climateLyrs[0].getXCor();
		lowerLeftY = climateLyrs[0].getYCor();
		cellSize = climateLyrs[0].getCellSize();
		
	}
	else if (geologyLyrCnt > 0)
	{
		totalRows = geologyLyrs[0].getNumOfRows();
		totalCols = geologyLyrs[0].getNumOfCols();
		lowerLeftX = geologyLyrs[0].getXCor();
		lowerLeftY = geologyLyrs[0].getYCor();
		cellSize = geologyLyrs[0].getCellSize();
	}
	else if (terrainLyrCnt > 0)
	{
		totalRows = terrainLyrs[0].getNumOfRows();
		totalCols = terrainLyrs[0].getNumOfCols();
		lowerLeftX = terrainLyrs[0].getXCor();
		lowerLeftY = terrainLyrs[0].getYCor();
		cellSize = terrainLyrs[0].getCellSize();
	}
	else if (vegeLyrCnt > 0)
	{
		totalRows = vegetationLyrs[0].getNumOfRows();
		totalCols = vegetationLyrs[0].getNumOfCols();
		lowerLeftX = vegetationLyrs[0].getXCor();
		lowerLeftY = vegetationLyrs[0].getYCor();
		cellSize = vegetationLyrs[0].getCellSize();
	}
	else if (otherLyrCnt > 0)
	{
		totalRows = otherLyrs[0].getNumOfRows();
		totalCols = otherLyrs[0].getNumOfCols();
		lowerLeftX = otherLyrs[0].getXCor();
		lowerLeftY = otherLyrs[0].getYCor();
		cellSize = otherLyrs[0].getCellSize();
	}
	  // obtain the x and y of the input sample
	/*
	图像中以左上角为起点，起点向右为x轴，逐渐增大，起点向下为y轴，逐渐减小。本程序中（见AscGrid.cpp）xCor，yCor为左上角的坐标。其中yCor在AscGrid.cpp中
	readAscGridGDAL等相关函数中被转换为了左下角的y坐标（xCor = pTransform[0];   yCor = pTransform[3] - pTransform[1]*totalRows;// 由左上角变为左下角）
	这时候我们的起点就是笛卡尔坐标系中的原点了，向右为x轴，逐渐增大，向上为y轴，逐渐增大。又左上角和左下角的x坐标是一样的，
	因此，xCor,yCor就可以当做是左下角坐标了即lowerLeftX，lowerLeftY
	我们在将样点的地理坐标转为图像中的行列位置时就以左下角为起点参考了。
	*/
    vector<string> sampleCoordinate;
    parseStr(string(sample), ',', sampleCoordinate);
    string sampleX = sampleCoordinate[0];
    string sampleY = sampleCoordinate[1];
    double curX = string_to_double(sampleX);
    double curY = string_to_double(sampleY);
    int rowIndex = totalRows - (int)((curY - lowerLeftY)/cellSize)-1;  // the position of Y in raster
    int colIndex = (int)((curX - lowerLeftX)/cellSize); // the position of X in raster

	sampleClimateV = new double[climateLyrCnt];
	sampleGeologyV = new double[geologyLyrCnt];
	sampleTerrainV = new double[terrainLyrCnt];
	sampleVegeV = new double[vegeLyrCnt];
	sampleOtherV = new double[otherLyrCnt];

    for (int i = 0; i < climateLyrCnt; i++)
    {
        sampleClimateV[i] = climateLyrs[i].values[rowIndex][colIndex];
    }
	
    for (int i = 0; i < geologyLyrCnt; i++)
    {
        sampleGeologyV[i] = geologyLyrs[i].values[rowIndex][colIndex];
    }
	for (int i = 0; i < terrainLyrCnt; i++)
	{
		sampleTerrainV[i] = terrainLyrs[i].values[rowIndex][colIndex];
	}
	for (int i = 0; i < vegeLyrCnt; i++)
	{
		sampleVegeV[i] = vegetationLyrs[i].values[rowIndex][colIndex];
	}
	for (int i = 0; i < otherLyrCnt; i++)
	{
		sampleOtherV[i] = otherLyrs[i].values[rowIndex][colIndex];
	}

	climateCnt = 0;
	geologyCnt = 0;
	terrainCnt = 0;
	vegeCnt = 0;
	otherCnt = 0;
	// load all layers using the function double *readAscGridGDAL_Block in AscGrid.cpp
	vector<string> environLyrs2;
	parseStr(string(environLyrsPath), '#', environLyrs2);
	for (int i = 0; i < typeOfLayers.size(); i++)
	{
		string type = typeOfLayers[i];
		AscGrid attributeLyr;
		double *rawData = attributeLyr.readAscGridGDAL_Block(environLyrs2[i], 0, 0, totalRows, totalCols);
		if (type == "Climate")
		{
			climate2DValues[climateCnt] = rawData;
			climateCnt++;
		}
		if (type == "Geology")
		{
			geology2DValues[geologyCnt] = rawData;
			geologyCnt++;
		}
		if (type == "Terrain")
		{
			terrain2DValues[terrainCnt] = rawData;
			terrainCnt++;
		}
		if (type == "Vegetation")
		{
			vege2DValues[vegeCnt] = rawData;
			vegeCnt++;
		}
		if (type == "Other")
		{
			other2DValues[otherCnt] = rawData;
		}
	} // end of for (int i = 0; i < typeOfLayers.size(); i++)

	double *climateSum;
	climateSum = new double[climateLyrCnt];
	int *climateCount;
	climateCount = new int[climateLyrCnt];
	double *geologySum;
	geologySum = new double[geologyLyrCnt];
	int *geologyCount;
	geologyCount = new int[geologyLyrCnt];
	double *terrainSum;
	terrainSum = new double[terrainLyrCnt];
	int *terrainCount;
	terrainCount = new int[terrainLyrCnt];
	double *vegeSum;
	vegeSum = new double[vegeLyrCnt];
	int *vegeCount;
	vegeCount = new int[vegeLyrCnt];
	double *otherSum;
	otherSum = new double[otherLyrCnt];
	int *otherCount;
	otherCount = new int[otherLyrCnt];

	climateSum = getSum(climate2DValues, climateLyrCnt, climateNoData, sampleClimateV, totalRows, totalCols);
	geologySum = getSum(geology2DValues, geologyLyrCnt, geologyNoData, sampleGeologyV, totalRows, totalCols);
	terrainSum = getSum(terrain2DValues, terrainLyrCnt, terrainNoData, sampleTerrainV, totalRows, totalCols);
	vegeSum = getSum(vege2DValues, vegeLyrCnt, vegeNoData, sampleVegeV, totalRows, totalCols);
	otherSum = getSum(other2DValues, otherLyrCnt, otherNoData, sampleOtherV, totalRows, totalCols);

	climateCount = getCount(climate2DValues, climateLyrCnt, climateNoData, totalRows, totalCols);
	geologyCount = getCount(geology2DValues, geologyLyrCnt, geologyNoData, totalRows, totalCols);
	terrainCount = getCount(terrain2DValues, terrainLyrCnt, terrainNoData, totalRows, totalCols);
	vegeCount = getCount(vege2DValues, vegeLyrCnt, vegeNoData, totalRows, totalCols);
	otherCount = getCount(other2DValues, otherLyrCnt, otherNoData, totalRows, totalCols);

	double **similarityValues;//相似度值
	double **uncertaintyValues;//不确定性值
	similarityValues = new double *[totalRows];
	uncertaintyValues = new double *[totalRows];
	double *similarityValuesStorage;
	double *uncertaintyValuesStorage;
	similarityValuesStorage = new double [totalRows*totalCols];
	uncertaintyValuesStorage = new double [totalRows*totalCols];
	for (int i = 0; i < totalRows; i++)
	{
		similarityValues[i] = &similarityValuesStorage[i*totalCols];
		uncertaintyValues[i] = &uncertaintyValuesStorage[i*totalCols];
	}

	ParallelSBMp inf(sampleClimateV, sampleGeologyV, sampleTerrainV,sampleVegeV, sampleOtherV,
		attributeRules,
		climateVRange,geologyVRange,terrainVRange,vegeVRange,otherVRange,
		rowIndex, colIndex, curX, curY, lowerLeftX, lowerLeftY, cellSize, totalRows,
		categoryIntegrationMethod, sampleIntegrationMethod,uncertaintyThreshold
		);
	inf.getPropertyMap(climateStd,geologyStd,terrainStd,vegeStd,otherStd,
		climate2DValues, geology2DValues,terrain2DValues, vege2DValues, other2DValues,
		similarityValues, uncertaintyValues,
		climateLyrCnt, geologyLyrCnt, terrainLyrCnt, vegeLyrCnt, otherLyrCnt,totalRows, totalCols,
		noData, climateNoData,geologyNoData,terrainNoData,vegeNoData,otherNoData,
		climateSum, geologySum, terrainSum, vegeSum, otherSum,
		climateCount, geologyCount, terrainCount, vegeCount, otherCount);//wtf

	NeighborSimi inf2(sampleClimateV, sampleGeologyV, sampleTerrainV,sampleVegeV, sampleOtherV,
		attributeRules,
		climateVRange,geologyVRange,terrainVRange,vegeVRange,otherVRange,
		rowIndex, colIndex, curX, curY, lowerLeftX, lowerLeftY, cellSize, totalRows,
		categoryIntegrationMethod, sampleIntegrationMethod,uncertaintyThreshold
		);
	inf2.getNeighborAvgSimi(climateStd,geologyStd,terrainStd,vegeStd,otherStd,
		climate2DValues, geology2DValues,terrain2DValues, vege2DValues, other2DValues,
		similarityValues, uncertaintyValues,
		climateLyrCnt, geologyLyrCnt, terrainLyrCnt, vegeLyrCnt, otherLyrCnt,totalRows, totalCols,
		noData, climateNoData,geologyNoData,terrainNoData,vegeNoData,otherNoData,
		climateSum, geologySum, terrainSum, vegeSum, otherSum,
		climateCount, geologyCount, terrainCount, vegeCount, otherCount);//wtf

	AscGrid similarityMap(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, similarityValues);
	AscGrid uncertaintyMap(totalCols, totalRows, lowerLeftX, lowerLeftY, cellSize, noData, uncertaintyValues);
	string similarityFile = simiLayerPath;
	string uncertaintyFile = uncertaintyLayerPath;

	/*similarityMap.createAscGridGADL(environLyrs[0],similarityFile);
	similarityMap.writeAscGridGDAL(similarityFile);
	uncertaintyMap.createAscGridGADL(environLyrs[0],uncertaintyFile);
	uncertaintyMap.writeAscGridGDAL(uncertaintyFile);*/


	delete[] similarityValues;
	delete[] uncertaintyValues;
	delete[] similarityValuesStorage;
	delete[] uncertaintyValuesStorage;
	delete[] climateNoData;
	delete[] geologyNoData;
	delete[] terrainNoData;
	delete[] vegeNoData;
	delete[] otherNoData;
	delete[] climateStd;
	delete[] geologyStd;
	delete[] terrainStd;
	delete[] vegeStd;
	delete[] otherStd;

	delete[] climateSum;
	delete[] climateCount;
	delete[] geologySum;
	delete[] geologyCount;
	delete[] terrainSum;
	delete[] terrainCount;
	delete[] vegeSum;
	delete[] vegeCount;
	delete[] otherSum;
	delete[] otherCount;

	// release memory

	if(climateLyrs != NULL)
		delete[] climateLyrs;
	if(geologyLyrs != NULL)
		delete[] geologyLyrs;
	if(terrainLyrs != NULL)
		delete[] terrainLyrs;
	if(vegetationLyrs != NULL)
		delete[] vegetationLyrs;
	if(otherLyrs != NULL)
		delete[] otherLyrs;

	if(climateVRange != NULL)
		delete[] climateVRange;
	if(geologyVRange != NULL)
		delete[] geologyVRange;
	if(terrainVRange != NULL)
		delete[] terrainVRange;
	if(vegeVRange != NULL)
		delete[] vegeVRange;
	if(otherVRange != NULL)
		delete[] otherVRange;

	if(sampleClimateV != NULL)
		delete[] sampleClimateV;
	if(sampleGeologyV != NULL)
		delete[] sampleGeologyV;
	if(sampleTerrainV != NULL)
		delete[] sampleTerrainV;
	if(sampleVegeV != NULL)
		delete[] sampleVegeV;
	if(sampleOtherV != NULL)
		delete[] sampleOtherV;

	if (climate2DValues != NULL)
	{
		delete[] climate2DValues;
	}
	if (geology2DValues != NULL)
	{
		delete[] geology2DValues;
	}
	if (terrain2DValues != NULL)
	{
		delete[] terrain2DValues;
	}
	if (vege2DValues != NULL)
	{
		delete[] vege2DValues;
	}
	if (other2DValues != NULL)
	{
		delete[] other2DValues;
	}
}