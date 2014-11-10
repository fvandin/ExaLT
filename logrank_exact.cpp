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
#include "logrank_exact.h"

using namespace std;

typedef struct Pair Pair;

struct Pair {
	double weight;
	double censor;
};

//global variables

//int n;

//#value fo the logrank statistic for which you want to compute the p-value
//global score 
double score; //not sure about the type yet

vector <double> weight_list;

//number of permutations that have logrank statistic >= score
int num_more_extreme;

//number of permutations that have logrank statistic < score
int num_less_extreme;

//functions
bool comp_pair(const Pair p1, const Pair p2) {
	return p1.weight>p2.weight;
}

bool comp_weight(const double w1, const double w2) {
	return w1>w2;
}

bool vec_contains(vector<int> v, int a) {
	int n_v = v.size();
	for (int i=0; i<n_v; i++) {
		if (v[i] == a)
			return true;
	}
	return false;
}


template <typename T>
bool next_combination(const T first, T k, const T last)
{
	/* Credits: Mark Nelson http://marknelson.us */
	if ((first == last) || (first == k) || (last == k))
		return false;
	T i1 = first;
	T i2 = last;
	++i1;
	if (last == i1)
		return false;
	i1 = last;
	--i1;
	i1 = k;
	--i2;
	while (first != i1)
	{
		if (*--i1 < *i2)
		{
			T j = k;
			while (!(*i1 < *j)) ++j;
			std::iter_swap(i1,j);
			++i1;
			++j;
			i2 = k;
			std::rotate(i1,j,last);
			while (last != j)
			{
				++j;
				++i2;
			}
			std::rotate(k,i2,last);
			return true;
		}
	}
	std::rotate(first,k,last);
	return false;
}

template <class T>
set<T> union_sets(set<T> a, set<T> b) {
	typename set<T>::iterator it = b.begin();
	while ( it != b.end() )  {
		a.insert(*it);
		it++;
	}
	return a;
}

template <class T>
vector<T> set2vector(set<T> s) {
	typename set<T>::iterator it = s.begin();
	vector<T> v;
	while ( it != s.end() ) {
		v.push_back(*it);
		it++;
	}
	return v;
}


template <class T>
vector<T> union_vectors(vector<T> a, vector<T> b) {
	typename vector<T>::iterator it = b.begin();
	while ( it != b.end() ) {
		a.push_back(*it);
		it++;
	}
	return a;
}


unsigned long long bcoeff(unsigned long long n, unsigned long long k) {
	if (k > n) {
		return 0;
	}
	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d) {
		r *= n--;
		r /= d;
	}
	return r;
}



double logrank(vector <int> input_list, int n1, vector <int> censoring_list) {
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
		/*
		#print "\ni: "+str(i)
		#print "observed: "+str(observed)
		#print "expected: "+str(expected)
		#print "curr. logrank_stat: "+str(logrank_stat)
		*/
	}
	return logrank_stat;

}




int exhaustive(vector<int> inlist, vector<int> censlist, double &p_value, double &p_value_left, double &p_value_right ){

/*	cout << "x= ";
	for (int i=0; i<inlist.size(); i++){
		cout << inlist[i];
	}
	cout << endl;
	cout << "c= ";
	for (int i=0; i<censlist.size(); i++){
                cout << censlist[i];
        }
        cout << endl;
*/
	int n_ones=0;

	for (int i=0; i<inlist.size(); i++) {
		if (inlist[i] == 1){ 
			n_ones += 1;
		}
	}

	double logstat = logrank(inlist,n_ones,censlist);

	//#build list to generate all possible combinations with same number of 0 and 1
	vector <int> pos_list;
	for (int i=0; i<inlist.size(); i++) {
		pos_list.push_back(i+1);
	}

       	unsigned long long num_more_extreme = 0;
	unsigned long long num_less_extreme = 0;
	unsigned long long more_extreme_right = 0;
	unsigned long long more_extreme_left = 0;

	vector<int> tmp_list;
	vector<int> pos_com_list;
	int pos_list_n = pos_list.size();
	do {
		pos_com_list.clear();
		for (int ii=0; ii<n_ones; ii++ ) {
			pos_com_list.push_back(*(pos_list.begin()+ii));
		}
		tmp_list.clear();		
		for (int i=0; i<pos_list_n; i++ ) {
				if ( vec_contains(pos_com_list, i+1)  )
					tmp_list.push_back(1);
				else
					tmp_list.push_back(0);
		}

		double tmp_logrank = logrank(tmp_list, n_ones, censlist);
		//print combination
		if (fabs(tmp_logrank) < fabs(logstat)) {
			num_less_extreme ++;
		}
		else {
			num_more_extreme ++;
			if (tmp_logrank >= fabs(logstat)){
				more_extreme_right++;
			}
			if (tmp_logrank <= -fabs(logstat)){
				more_extreme_left++;
			}
		}

	} while( next_combination( pos_list.begin(), pos_list.begin()+n_ones, pos_list.end() ) );

	p_value = double(num_more_extreme)/double(num_more_extreme+num_less_extreme);
	p_value_left = double(more_extreme_left)/double(num_more_extreme+num_less_extreme);
	p_value_right = double(more_extreme_right)/double(num_more_extreme+num_less_extreme);
}

