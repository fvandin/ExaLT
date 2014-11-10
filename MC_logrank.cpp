/*Copyright 2013 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose other than its incorporation into a
commercial product is hereby granted without fee, provided that the
above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

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

double logrank_MC(vector <int> input_list, int n1, vector <int> censoring_list) {
        double curr_population = (double)input_list.size();
        double curr_ones = double(n1);
        double logrank_stat = 0.0;
        int n = input_list.size();
        double observed, expected;
        for (int i=0; i<n; i++) {
                observed = (double)input_list[i];
	        expected = curr_ones/curr_population;
	    	if (input_list[i] == 1)
                       curr_ones -= 1.0;
		curr_population -= 1.0;
		logrank_stat += (observed - expected)*censoring_list[i];
	}
	return logrank_stat;
}

int MC_estimate( vector<int> inlist, vector<int> censlist, double &p_value, double &p_value_left, double & p_value_right, unsigned long long num_iter ) {

	int n_ones = 0;
        for (int i=0; i<inlist.size(); i++) {
                if (inlist[i] == 1) {
                        n_ones += 1;
		}
	}

	double logstat = logrank_MC(inlist,n_ones,censlist);

	//initialize random seed
	srand( unsigned ( time(0) ) );

	//sample random permutations at random, and use that to computer p-val
	unsigned long long found = 0;
	unsigned long long found_left = 0;
	unsigned long long found_right = 0;
	for (int i =1; i<= num_iter; i++){
		//get random permutation
		random_shuffle ( inlist.begin(), inlist.end() );
		//get logrank statistic of current permutation
		double curr_logstat = logrank_MC(inlist, n_ones,censlist);
		if (fabs( curr_logstat ) >= fabs(logstat))
			found++;
		if (curr_logstat >= fabs(logstat)){
			found_right++;
		}
		if (curr_logstat <= -fabs(logstat)){
			found_left++;
		}
	}
	//estimated p-val
	double p_val;
	p_value = (double)found/((double) num_iter);
	p_value_left = (double)found_left/((double) num_iter);
	p_value_right = (double)found_right/((double) num_iter);
}
