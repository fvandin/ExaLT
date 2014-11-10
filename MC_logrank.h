#include <limits>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <utility>


#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using namespace std;

int MC_estimate( vector<int> inlist, vector<int> censlist, double &p_value, double &p_value_left, double &p_value_right, unsigned long long num_iter );
