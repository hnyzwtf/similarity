#include <sstream>
#include "NeighborSimi.h"
#include "Distance.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
using namespace std;

double NeighborSimi::string_to_double( const std::string& s )
{
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 }

NeighborSimi::NeighborSimi(){

}

NeighborSimi::NeighborSimi(double * sampleClimateV, double * sampleGeologyV, double * sampleTerrainV,double * sampleVegetationV, double * sampleOtherV,
			vector<string> &AttriRules,
             double * climateVRange,double * geologyVRange,double * terrainVRange,double * vegeVRange,double * otherVRange,
			 int rowIndex, int colIndex, double curX, double curY, double lowerLeftX, double lowerLeftY, double cellSize, int totalRows,
             string catIntegrationMethod, string sampleIntegrationMethod,  string uncertaintyThreshold
			 ):attriRules(AttriRules.size()),catIntegration(catIntegrationMethod),sampleIntegration(sampleIntegrationMethod)
{
     this->sampleClimateV = sampleClimateV;
     this->sampleGeologyV = sampleGeologyV;
     this->sampleTerrainV = sampleTerrainV;
	 this->sampleVegetationV = sampleVegetationV;
	 this->sampleOtherV = sampleOtherV;

	 this->climateVRange = climateVRange;
	 this->geologyVRange = geologyVRange;
	 this->terrainVRange = terrainVRange;
	 this->vegeVRange = vegeVRange;
	 this->otherVRange = otherVRange;

	 this->rowIndex = rowIndex;
	 this->colIndex = colIndex;
	 this->curX = curX;
	 this->curY = curY;
	 this->lowerLeftX = lowerLeftX;
	 this->lowerLeftY = lowerLeftY;
	 this->cellSize = cellSize;
	 this->totalRows = totalRows;

	 for (int i = 0; i < AttriRules.size(); i ++){
			vector<string> ruleParameterArray;
			split(AttriRules[i],"?",ruleParameterArray);
			this->attriRules[i].setCategory(ruleParameterArray[0]);
			this->attriRules[i].setDistCalculationMethod(ruleParameterArray[1]);
	}
	this->uncertaintyThreshold = atof(uncertaintyThreshold.c_str());
}

void NeighborSimi::getNeighborAvgSimi(double *climateStd,double *geoStd,double *terrainStd,double *vegetationStd,double *otherStd,
					double ** climate2Dvalues, double ** geo2Dvalues,double ** terrain2Dvalues, double ** vege2Dvalues, double ** other2Dvalues,
				   double ** similarityVals, double ** uncertaintyVals,
				   int ClimateLyrCnt, int GeogLyrCnt, int TerrainLyrCnt, int VegeLyrCnt, int OtherLyrCnt,int rows, int cols, 
				   double noData,double *climateNoData, double *geoNoData,double *terrainNoData, double *vegetationNoData,double *otherNoData,
				  double *climateSum, double *geologySum, double *terrainSum, double *vegetationSum, double *otherSum,
				int *climateCount, int *geologyCount, int *terrainCount, int *vegetationCount, int *otherCount)//wtf
   {
	   AttributeRule* climateRules = new AttributeRule[ClimateLyrCnt];
		AttributeRule* geologyRules = new AttributeRule[GeogLyrCnt];
		AttributeRule* terrainRules = new AttributeRule[TerrainLyrCnt];
		AttributeRule* vegetationRules = new AttributeRule[VegeLyrCnt];
		AttributeRule* otherRules = new AttributeRule[OtherLyrCnt];
		int  numOfLyrs = ClimateLyrCnt + GeogLyrCnt + TerrainLyrCnt + VegeLyrCnt + OtherLyrCnt;
		int climateCnt =0;
		int geologyCnt = 0;
		int terrainCnt =0;
		int vegeCnt =0;
		int otherCnt =0;       
		for(int attriDataLyrIdx = 0; attriDataLyrIdx <numOfLyrs; attriDataLyrIdx++){

			string category = attriRules[attriDataLyrIdx].getCategory();

			if(category == "Climate"){
				 climateRules[climateCnt] = attriRules[attriDataLyrIdx];
				 climateCnt ++;
			}else if(category == "Geology"){
				  geologyRules[geologyCnt] = attriRules[attriDataLyrIdx];
				  geologyCnt ++;
			}else if(category == "Terrain"){
				  terrainRules[terrainCnt] = attriRules[attriDataLyrIdx];
				  terrainCnt ++;
			}else if(category == "Vegetation"){
				  vegetationRules[vegeCnt] = attriRules[attriDataLyrIdx];
				  vegeCnt ++;

			}else if(category == "Other"){
				   otherRules[otherCnt] = attriRules[attriDataLyrIdx];
				   otherCnt ++;
			}

		}

		double* CurrentCellAttributeSimilarities = new double[5];
		int* NumOfLayersInCategories = new int[5];
		double* ClimateSimilarities = new double[ClimateLyrCnt];
		double* GeologySimilarities = new double[GeogLyrCnt];
		double* TerrainSimilarities = new double[TerrainLyrCnt];
		double* VegetationSimilarities = new double[VegeLyrCnt];
		double* OthersSimilarities = new double[OtherLyrCnt];	

		double simiPointsX;// �������ɵĿ�����������
		double simiPointsY;
		Distance distance;
		int sampleNeighborCnt = 0; //��������������Χһ�������ڵ���Ԫ����
		double sampleNeighborSum = 0; //��������������Χһ��������������Ԫ���ƶ����
		double sampleNeighborAvg = 0; //��������������Χx����������Ԫ���ƶ����ֵ
		int candidateNeighborCnt = 0; //�����ѡ���㼴�ɹ�ѡ���������Χ�����ڵ���Ԫ����		
		double candidateSum = 0; //ÿ����ѡ������Χ������Ԫ���ƶȵĺ�	
		vector <double> candidateNeighborAvg;//����ÿ����ѡ������Χ������Ԫ�����ƶ���ƽ������ÿ����ѡ�������ƶ�ƽ��ֵ��vector
		vector<int> candidateRow;//��ѡ������ͼ���е�row�����к�
		vector<int> candidateCol;
//=================================================================================================================================
		for(int rowIdx = 0; rowIdx < rows; rowIdx++)
		{
			 for(int colIdx = 0; colIdx < cols; colIdx++)
			 {
				
					double MaxSimilarity = 0;
					// climate
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < ClimateLyrCnt; attriDataLyrIdx++){
						double gridV = climate2Dvalues[attriDataLyrIdx][rowIdx * cols + colIdx];
						if(abs(gridV - climateNoData[attriDataLyrIdx])< VERY_SMALL){
							ClimateSimilarities[attriDataLyrIdx] = -9999;
						}else{
								double sampleV = sampleClimateV[attriDataLyrIdx];//դ��ͼ��������λ���ϵ���ֵ
								double range = climateVRange[attriDataLyrIdx];
								ClimateSimilarities[attriDataLyrIdx] = climateRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,climateSum[attriDataLyrIdx],climateCount[attriDataLyrIdx],climateStd[attriDataLyrIdx]);                                 
							}
					}
					double FinalClimateSimilarity = catIntegration.GetCategorySimilarity(ClimateSimilarities, ClimateLyrCnt);
					CurrentCellAttributeSimilarities[0] = FinalClimateSimilarity;
					NumOfLayersInCategories[0] = ClimateLyrCnt;
					// geology
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < GeogLyrCnt; attriDataLyrIdx++){                      
						double gridV = geo2Dvalues[attriDataLyrIdx][rowIdx * cols + colIdx];
						if(abs(gridV - geoNoData[attriDataLyrIdx])< VERY_SMALL){
							GeologySimilarities[attriDataLyrIdx] = -9999;
						}else{

							double sampleV = sampleGeologyV[attriDataLyrIdx];
							double range = geologyVRange[attriDataLyrIdx];
							GeologySimilarities[attriDataLyrIdx] = geologyRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,geologySum[attriDataLyrIdx],geologyCount[attriDataLyrIdx],geoStd[attriDataLyrIdx]);
						}

					}

					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(GeologySimilarities, GeogLyrCnt);
					CurrentCellAttributeSimilarities[1] = FinalClimateSimilarity;
					NumOfLayersInCategories[1] = GeogLyrCnt;
					//terrain
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < TerrainLyrCnt; attriDataLyrIdx++){
						double gridV = terrain2Dvalues[attriDataLyrIdx][rowIdx * cols + colIdx];
						if(abs(gridV - terrainNoData[attriDataLyrIdx])< VERY_SMALL){
							TerrainSimilarities[attriDataLyrIdx] = -9999;							
						}else{
							double sampleV = sampleTerrainV[attriDataLyrIdx];
							double range = terrainVRange[attriDataLyrIdx];
							TerrainSimilarities[attriDataLyrIdx] = terrainRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,terrainSum[attriDataLyrIdx],terrainCount[attriDataLyrIdx],terrainStd[attriDataLyrIdx]);
						}
					}

					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(TerrainSimilarities, TerrainLyrCnt);                    
					CurrentCellAttributeSimilarities[2] = FinalClimateSimilarity;
					NumOfLayersInCategories[2] = TerrainLyrCnt;
					//vegetation
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < VegeLyrCnt; attriDataLyrIdx++){                 
						double gridV = vege2Dvalues[attriDataLyrIdx][rowIdx * cols + colIdx];
						if(abs(gridV - vegetationNoData[attriDataLyrIdx])< VERY_SMALL){
							VegetationSimilarities[attriDataLyrIdx] = -9999;
						}else{
							double sampleV = sampleVegetationV[attriDataLyrIdx];
							double range = vegeVRange[attriDataLyrIdx];							
							VegetationSimilarities[attriDataLyrIdx] = vegetationRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,vegetationSum[attriDataLyrIdx],vegetationCount[attriDataLyrIdx],vegetationStd[attriDataLyrIdx]);
						}

					}

					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(VegetationSimilarities, VegeLyrCnt);

					CurrentCellAttributeSimilarities[3] = FinalClimateSimilarity;
					NumOfLayersInCategories[3] = VegeLyrCnt;
					// others
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < OtherLyrCnt; attriDataLyrIdx++){                        
						double gridV = other2Dvalues[attriDataLyrIdx][rowIdx * cols + colIdx];
						if(abs(gridV - otherNoData[attriDataLyrIdx])< VERY_SMALL){
							OthersSimilarities[attriDataLyrIdx] = -9999;
						}else{

							double sampleV = sampleOtherV[attriDataLyrIdx];
							double range = otherVRange[attriDataLyrIdx];							
							OthersSimilarities[attriDataLyrIdx] = otherRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,otherSum[attriDataLyrIdx],otherCount[attriDataLyrIdx],otherStd[attriDataLyrIdx]);
						}

					}
					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(OthersSimilarities, OtherLyrCnt);
					CurrentCellAttributeSimilarities[4] = FinalClimateSimilarity;
					NumOfLayersInCategories[4] = OtherLyrCnt;

					int count = 5;
					//�˺��������ۺ�ĳ��������ĳ�����Ʋ��Ķ���������ƶ�
					double similarity = sampleIntegration.GetSampleSimilarity(CurrentCellAttributeSimilarities, count, NumOfLayersInCategories);                    			
				
					if (similarity >= MaxSimilarity)
					{
						MaxSimilarity = similarity;
					}							

				double uncertainty = 1 - MaxSimilarity;
				if ( uncertainty < this->uncertaintyThreshold ){
					uncertaintyVals[rowIdx][colIdx] = uncertainty;
					similarityVals[rowIdx][colIdx] = similarity;
			
				}else{
					uncertaintyVals[rowIdx][colIdx] = noData;
					similarityVals[rowIdx][colIdx] = noData;
				}
				// ��������������Χһ������Χ�����ƶ�ֵ��ƽ��
				// ����б�ʾ��������������Χ���������к����Ϊ8��դ����Ԫ�������ֵʱ���ܰ���������������ƶ�ֵ����Χ�����ƶ�ֵ��Ϊ��
				if (abs(rowIdx-rowIndex) < 8 && abs(colIdx-colIndex) < 8 && (rowIdx != rowIndex) && (colIdx != colIndex)
					&& similarityVals[rowIdx][colIdx] != noData)
				{			
					sampleNeighborCnt++;				
					sampleNeighborSum += similarityVals[rowIdx][colIdx];
				}
				// rowIndex��colIndex������������λ�ã�rowIdx��colIdx��դ��ͼ����ÿһ��դ�������к�
				// ���ͼ��������һ������ƶ�ֵ��������������ƶ�ֵ֮��С��0.02���ͰѴ˵���Ϊ��ѡ����
				if (abs(similarityVals[rowIndex][colIndex] - similarityVals[rowIdx][colIdx]) < 0.02)
				{
					simiPointsX = lowerLeftX + (double)(colIdx * cellSize);//��ѡ��������Ӧ������x
					simiPointsY = lowerLeftY + (double)(totalRows - rowIdx - 1) * cellSize;
					candidateRow.push_back(rowIdx);//�����к�ѡ������ͼ���ϵ����кŴ���vector������
					candidateCol.push_back(colIdx);
					//cout<<"These are similar points: "<<rowIdx<<","<<colIdx<<endl;	
				}						
			}//end of for(int colIdx = 0; colIdx < cols; colIdx++)
		}// end of for(int rowIdx = 0; rowIdx < rows; rowIdx++)
//=================================================================================================================================
		sampleNeighborAvg = sampleNeighborSum/sampleNeighborCnt;
		
		cout<<"this is sampleNeighborAvg   "<<sampleNeighborAvg<<endl;

		for (int i = 0; i < candidateRow.size(); i++)// candidateRow.size()���Ǻ�ѡ����ĸ���
		{	
			candidateNeighborCnt = 0;
			candidateSum = 0;
			for(int rowIdx = 0; rowIdx < rows; rowIdx++)
			{
				for(int colIdx = 0; colIdx < cols; colIdx++)
				{
					if (abs(rowIdx-candidateRow[i]) < 8 && abs(colIdx-candidateCol[i]) < 8 && (rowIdx != candidateRow[i]) && (colIdx != candidateCol[i])
						&& similarityVals[rowIdx][colIdx] != noData)
					{			
						candidateNeighborCnt++;				
						candidateSum += similarityVals[rowIdx][colIdx];
					}
				}
			}//candidateNeighborAvg��ŵľ���ÿ����ѡ������Χդ�����ƶȵ�ƽ��ֵ���ж��ٸ���ѡ����,���ж��ٸ���ֵ			
			candidateNeighborAvg.push_back(candidateSum / candidateNeighborCnt);
		}	
		for (int i = 0; i < candidateNeighborAvg.size(); i++)
		{
			cout<<"this is candidateNeighborAvg--->"<<candidateNeighborAvg[i]<<endl;
		}

		delete[] climateRules;
		delete[] geologyRules;
		delete[] terrainRules;
		delete[] vegetationRules;
		delete[] otherRules;
		delete[] ClimateSimilarities;
		delete[] GeologySimilarities;
		delete[] TerrainSimilarities;
		delete[] VegetationSimilarities;
		delete[] OthersSimilarities;
		delete [] CurrentCellAttributeSimilarities;
		delete [] NumOfLayersInCategories;
   }

void NeighborSimi::split(string s, string sep, vector<string> &flds){
	if (s.length() == 0)
		return;
	string fld;
	unsigned int i = 0, j = 0;
	do {
		j = s.find_first_of(sep,i);
		if (j > s.length())
			j = s.length();
		//fld = string(s,i,j-i);
		fld = s.substr(i,j-i);
		flds.push_back(fld);
		i = j + 1;
	}while(j < s.length());
}

