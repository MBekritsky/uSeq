#include "HelperFunc.h"

using namespace std;

vector< bool >  Dec2Bin (int dec)
{
	vector< bool > bin;
	
	for(unsigned int i = 0; dec >= 1; i++)
		{
			bin . push_back(dec % 2);
			dec /= 2;
		}
	
	reverse(bin . begin(), bin . end());	
	
	return bin;
}

unsigned int Bin2Dec (vector< bool > bin)
{
	unsigned int num = 0;
	
	vector< bool > tbin = bin;
	
	reverse(tbin . begin(), tbin . end());
	
	for(int i = 0; i < (int) tbin . size(); i++)
		{
			num += (tbin[i] * ((int)  pow(2.0, i)));
		}
	return num;

}
