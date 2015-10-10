#include "Distance.h"
using namespace std;

Distance::Distance()
{

}
double Distance::getEucliDistance(double sampleX, double sampleY, double simiPointsX, double simiPointsY)
{
	double d = sqrt(pow((sampleX-simiPointsX),2) + pow((sampleY-simiPointsY),2));
	return d;
}