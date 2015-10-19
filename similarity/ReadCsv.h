#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
class ReadCsv
{
public:
	ReadCsv();
	string Trim(string& str);
	void getCsvContent(string csvPath);
	double string_to_double( const std::string& s );
	vector <double> getCsvX();
	vector <double> getCsvY();
protected:
private:
	vector<double> csvX;
	vector<double> csvY;	
};