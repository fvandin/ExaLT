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

#ifndef TABLE_H
#define TABLE_H

#include <sstream>
#include <fstream>
#include <string>
#include <map>

#include <iostream>
#include <vector>
#include <algorithm>

#include <stdlib.h>

using namespace std;

#define LARGE_TO_SMALL 0
#define SMALL_TO_LARGE 1


struct order_pair {
	double entry;
	int index;
};

/*
bool compareS2B(struct order_pair a, struct order_pair b) {
	return a.entry<b.entry;
}
bool compareB2S(struct order_pair a, struct order_pair b) {
	return a.entry>b.entry;
}
*/

int parseHeader(string line, vector<string> *out);

bool compareS2B(struct order_pair a, struct order_pair b);
bool compareB2S(struct order_pair a, struct order_pair b);

class table {
	private:
		vector< vector<string> > data;
		map<string,int> colLookUpTable;
		vector<string> title;
	public:
		table(char* filename);
		void addTitle( vector<string> title_row);
		int addRow(vector<string> row);
		vector<string> getACol(int i);
		int rows();
		int cols();
		vector<int> getAColAsInt(int i);
		vector<int> vec_str2int(vector<string> v);
		//sort the table base of ith column, wayToSort:LARGE_TO_SMALL or SMALL_TO_LARGE
		//output an order sequence
		vector<int> sort_order(int j, int wayToSort); 
		vector<int> getColAsIntByOrder(int j, vector<int> order);
		vector<int> getColAsIntByOrder(string name, vector<int> order);
		vector<int> getStatusByOrder(int j, vector<int> order);
		vector<int> getStatusByOrder(string st, vector<int> order);
};

//vector utilities
void disp_vec_int(vector<int> v, string name);

void disp_vec_string(vector<string> v, string name);


#endif
