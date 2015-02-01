/*
HelperFunc.h

ABSTRACT: random functions that don't fit into any other categories

CREATION DATE: 03.05.2011

LAST REVISION: 03.05.2011

AUTHOR:
Mitchell Bekritsky
*/

#ifndef __HELPER_FUNC_H__
#define __HELPER_FUNC_H__

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

vector< bool > Dec2Bin (int dec);
unsigned int Bin2Dec (vector< bool > bin);

#endif /* HelperFunc.h */
