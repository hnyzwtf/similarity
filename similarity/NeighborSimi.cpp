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

		double simiPointsX;// 程序生成的可替代点的坐标
		double simiPointsY;
		Distance distance;
		
		int candidateNeighborCnt = 0; //距离候选样点即可供选择的样点周围邻域内的像元个数		
		double candidateSum = 0; //每个候选样点周围所有像元相似度的和	
		vector <double> candidateNeighborAvg;//距离每个候选样点周围所有像元的相似度求平均，是每个候选样点相似度平均值的vector
		vector<int> candidateRow;//候选样点在图像中的row，即行号
		vector<int> candidateCol;
		vector <double> candidateVariance;//候选样点的方差集合
		double candidateVarianceCnt = 0;
		double candidateVarianceSum = 0;
		
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
								double sampleV = sampleClimateV[attriDataLyrIdx];//栅格图像在样点位置上的数值
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
					//此函数用来综合某个样点与某个待推测点的多个类别的相似度
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
						
				// rowIndex和colIndex是输入的样点的位置，rowIdx和colIdx是栅格图像上每一个栅格点的行列号
				// 如果图像中任意一点的相似度值和输入样点的相似度值之差小于0.02，就把此点作为候选样点
				if (abs(similarityVals[rowIndex][colIndex] - similarityVals[rowIdx][colIdx]) < 0.02)
				{
					simiPointsX = lowerLeftX + (double)(colIdx * cellSize);//候选样点所对应的坐标x
					simiPointsY = lowerLeftY + (double)(totalRows - rowIdx - 1) * cellSize;
					candidateRow.push_back(rowIdx);//将所有候选样点在图像上的行列号存入vector容器中
					candidateCol.push_back(colIdx);
					//cout<<simiPointsX<<","<<simiPointsY<<endl;	
				}						
			}//end of for(int colIdx = 0; colIdx < cols; colIdx++)
		}// end of for(int rowIdx = 0; rowIdx < rows; rowIdx++)
//=================================================================================================================================
	
		for (int i = 0; i < candidateRow.size(); i++)// candidateRow.size()就是候选样点的个数
		{	
			candidateNeighborCnt = 0;
			candidateSum = 0;
			for(int rowIdx = 0; rowIdx < rows; rowIdx++)
			{
				for(int colIdx = 0; colIdx < cols; colIdx++)
				{
					// 语句中表示的是候选样点周围和此候选行列号相差为8的栅格像元，计算时不能包含候选样点的相似度值
					if (abs(rowIdx-candidateRow[i]) < 10 && abs(colIdx-candidateCol[i]) < 10 && (rowIdx != candidateRow[i]) && (colIdx != candidateCol[i])
						&& similarityVals[rowIdx][colIdx] != noData)
					{			
						candidateNeighborCnt++;				
						candidateSum += similarityVals[rowIdx][colIdx];
					}
				}
			}//candidateNeighborAvg存放的就是每个候选样点周围栅格相似度的平均值，有多少个候选样点,就有多少个均值			
			candidateNeighborAvg.push_back(candidateSum / candidateNeighborCnt);
		}	
		for (int i = 0; i < candidateNeighborAvg.size(); i++)
		{
			cout<<"候选点周围邻域像元相似度的平均--->"<<candidateNeighborAvg[i]<<endl;
		}
		// 计算每个候选样点周围邻域相似度的方差
		for (int i = 0; i < candidateRow.size(); i++)
		{	
			candidateVarianceCnt = 0;
			candidateVarianceSum = 0;
			for(int rowIdx = 0; rowIdx < rows; rowIdx++)
			{
				for(int colIdx = 0; colIdx < cols; colIdx++)
				{
					if (abs(rowIdx-candidateRow[i]) < 10 && abs(colIdx-candidateCol[i]) < 10 && (rowIdx != candidateRow[i]) && (colIdx != candidateCol[i])
						&& similarityVals[rowIdx][colIdx] != noData)
					{			
						candidateVarianceCnt++;				
						candidateVarianceSum += pow((similarityVals[rowIdx][colIdx] - candidateNeighborAvg[i]),2) ;
					}
				}
			}//candidateNeighborAvg存放的就是每个候选样点周围栅格相似度的平均值，有多少个候选样点,就有多少个均值			
			candidateVariance.push_back(candidateVarianceSum / candidateVarianceCnt);
		}	
		cout<<"-----------------------------------------------------"<<endl;
		for (int i = 0; i < candidateVariance.size(); i++)
		{
			cout<<"候选点周围邻域像元相似度的方差--->"<<candidateVariance[i]<<endl;
		}
		// 对均值和方差结果进行排序，其中均值按照从大到小的顺序，如果有几个元素均值相等的话，他们对应的方差按从小到大排序
		int candidateNum = candidateRow.size();
		double temp;
		for (int i = 1; i < candidateNum; i++)
		{
			for (int j = 0; j < candidateNum; j++)
			{
				if (candidateNeighborAvg[i] > candidateNeighborAvg[j])
				{
					temp = candidateNeighborAvg[i];
					candidateNeighborAvg[i] = candidateNeighborAvg[j];
					candidateNeighborAvg[j] = temp;
					// 因为candidateRow,candidateCol,avg,variance里面的元素都要对应，因此，排序时平均值每个元素的位置发生了变化，相应的row,col,variance
					//位置也要发生变化
					temp = candidateVariance[i];
					candidateVariance[i] = candidateVariance[j];
					candidateVariance[j] = temp;

					temp = candidateRow[i];
					candidateRow[i] = candidateRow[j];
					candidateRow[j] = temp;

					temp = candidateCol[i];
					candidateCol[i] = candidateCol[j];
					candidateCol[j] = temp;
				}
				if (abs(candidateNeighborAvg[i] - candidateNeighborAvg[j]) < 0.00001)//如果两个元素相同
				{
					if (candidateVariance[i] < candidateVariance[j])//对方差进行从小到大的排序
					{
						temp = candidateNeighborAvg[i];
						candidateNeighborAvg[i] = candidateNeighborAvg[j];
						candidateNeighborAvg[j] = temp;

						temp = candidateVariance[i];
						candidateVariance[i] = candidateVariance[j];
						candidateVariance[j] = temp;

						temp = candidateRow[i];
						candidateRow[i] = candidateRow[j];
						candidateRow[j] = temp;

						temp = candidateCol[i];
						candidateCol[i] = candidateCol[j];
						candidateCol[j] = temp;
					}
				}
			}

		}
		for (int i = 0; i < candidateNum; i++)
		{
			cout<<candidateNeighborAvg[i]<<",";
		}
		cout<<endl;
		cout<<"------------------------------------------------------------------------------------"<<endl;
		for (int i = 0; i < candidateNum; i++)
		{		
			cout<<candidateVariance[i]<<",";
		}
		cout<<"------------------------------------------------------------------------------------"<<endl;
		for (int i = 0; i < candidateNum; i++)
		{
			cout<<"("<<candidateRow[i]<<","<<candidateCol[i]<<")"<<",";
		}
	/*	int index = 0;
		double min = 9999;
		for (int i = 0; i < candidateVariance.size(); i++)
		{		
			if (abs(sampleVariance - candidateVariance[i]) < min)
			{
				min = abs(sampleVariance - candidateVariance[i]);
				index = i;
			}
		}*/
		

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


