#include "ReadCsv.h"
using namespace std;

ReadCsv::ReadCsv()
{

}
double ReadCsv::string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}
string ReadCsv::Trim(string& str)
{
	str.erase(0,str.find_first_not_of(" \t\r\n"));
	str.erase(str.find_last_not_of(" \t\r\n") + 1);
	return str;
}
void ReadCsv::getCsvContent(string csvPath)
{
	ifstream fin(csvPath);
	string line;    
	while (getline(fin, line)) {
		//cout << line << endl;
		istringstream sin(line);    
		vector<string> fields;    
		string field;
		while (getline(sin, field, ',')) {
			fields.push_back(field);    
		}
		//string name = Trim(fields[0]);  
		string x = Trim(fields[1]);   
		string y = Trim(fields[2]);  
		if ((x != "X") && (y != "Y"))
		{
			double xx = string_to_double(x);
			double yy = string_to_double(y);
			csvX.push_back(xx);
			csvY.push_back(yy);
		}			
	}
}
vector <double> ReadCsv::getCsvX()
{
	return csvX;
}
vector <double> ReadCsv::getCsvY()
{
	return csvY;
}