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

#include "table.h"


table::table(char* filename) {
	string line;
	ifstream infile(filename);
	//table t;

	bool title_built = false;
	//build table
	while (getline(infile, line))  // this does the checking!
	{
		//cout<<line<<endl<<endl<<endl;
		vector<string> out;
		parseHeader(line, &out);
		if (title_built == false) {
			addTitle(out);
			title_built=true;
#ifdef DBG
			cout<<"add title"<<endl;
#endif
		}
		else {
#ifdef DBG
			disp_vec_string(out, "out");
#endif
			addRow(out);
#ifdef DBG
			cout<<"add row"<<endl;
#endif
		}
#ifdef DBG
		//cout<<"!!!!!!!!!! tuples: "<<t.rows()<<endl;
		cout<<endl<<endl;
#endif
	}


#ifdef DBG
	cout<<"In table:table,  list the whole table"<<endl;
	map<string,int>::iterator it2;  
	for (it2= colLookUpTable.begin(); it2!=colLookUpTable.end(); it2++) {
		cout<<it2->first<<", "<<it2->second<<endl;
	}
#endif
}


void table::addTitle( vector<string> title_row) {
	int n = title_row.size();
	for (int i=0; i<n; i++) {
		title.push_back(title_row[i]);
		//build col lookup table as well
		colLookUpTable[title_row[i]] = i;
#ifdef DBG
		cout<<"In table::addTitle, colLookUpTable["<<title_row[i]<<"] = "<<i<<endl;
#endif
	}
	//test code
	//cout<<"colLookUpTable[\"WHAT\"] = "<<colLookUpTable["WHAT"]<<endl;
	/*
	map<string, int>::iterator it = colLookUpTable.find("what");
	if (it==colLookUpTable.end()) {
		cout<<"no found"<<endl;
	} else {
		cout<<"it->first="<<it->first<<endl;
		cout<<"it->second="<<it->second<<endl;
	}*/
}

int table::addRow(vector<string> row) {
	//check number of attributes
	if (row.size() != cols() ) {
		cerr<<"Error in table::addRow: The row to be inserted, "<<endl;
		disp_vec_string(row, "row");
		cerr<<"has "<<row.size()<<" columns, but the schema has "<<cols()<<" columns."<<endl;
		return -1;
	}
	data.push_back(row);
#ifdef DBG
	cout<<"in table::addRow data.size()="<<data.size()<<endl;
#endif
	return 0;
}

vector<string> table::getACol(int i) {
	//error checking
	vector<string> ret;
	int r = rows();
	for (int j=0; j<r; j++) {
		ret.push_back(data[j][i]);
	}
	return ret;
}

int table::rows() {
	return data.size();
}

int table::cols() {
	//if ( rows() > 0 )
	return title.size();
}

vector<int> table::getAColAsInt(int i) {
	return vec_str2int( getACol(i) );
}

vector<int> table::vec_str2int(vector<string> v) {
	vector<int> ret;
	int n=v.size();
	for(int i=0; i<n; i++) {
		ret.push_back( atoi(v[i].c_str()) );
	}
	return ret;
}

vector<int> table::sort_order(int j, int wayToSort) {
	//build the order pair
	int r = rows();
#ifdef DBG
	cout<<"in table::sort_order, rows="<<r<<endl;
#endif
	vector<order_pair> order_vec;
	for (int i=0; i<r; i++) { //for all rows
		struct order_pair tmp;
#ifdef DBG
		cout<<"in table::sort_order, data["<<i<<"]["<<j<<"]="<<data[i][j]<<endl;
#endif
		tmp.entry = atof(data[i][j].c_str());
		tmp.index=i;
		order_vec.push_back(tmp);
	}

	//sort
	if (wayToSort==SMALL_TO_LARGE)
		sort(order_vec.begin(), order_vec.end(), compareS2B);
	else
		sort(order_vec.begin(), order_vec.end(), compareB2S);

	//return the order vec
	vector<int> index_vec;
	for (int i=0; i<r; i++) { //for all rows
		index_vec.push_back(order_vec[i].index);
	}
	return index_vec;
}


vector<int> table::getColAsIntByOrder(int j, vector<int> order) {
	int r = rows();
#ifdef DBG
	cout<<"in table::getColAsIntByOrder, r="<<r<<endl;
	cout<<"in table::getColAsIntByOrder, cols="<<cols()<<endl;
	cout<<"in table::getColAsIntByOrder, j="<<j<<endl;

	cout<<"data.size()="<<data.size()<<endl;
	cout<<"data[0].size()="<<data[0].size()<<endl;
	cout<<"data[1].size()="<<data[1].size()<<endl;

	cout<<"data[0][4]="<<data[0][4]<<endl;
	cout<<"data[1][4]="<<data[1][4]<<endl;
#endif
	vector<int> ret;
	for (int i=0; i<r; i++) {
#ifdef DBG
		cout<<"in table::getColAsIntByOrder, i="<<i<<endl;
		cout<<"in table::getColAsIntByOrder, order["<<i<<"]="<<order[i]<<endl;
#endif
		ret.push_back( atoi(data[order[i]][j].c_str()) );
	}
	return ret;
}


vector<int> table::getColAsIntByOrder(string name, vector<int> order) {
	//look up index of name
	//cout<<"in getColAsIntByOrder(name version), colLookUpTable["<<name<<"]="<<colLookUpTable[name]<<endl;
	map<string,int>::iterator it = colLookUpTable.find(name);

#ifdef DBG
	map<string,int>::iterator it2;  
	cout<<"in table::getColAsIntByOrder,  list the whole table"<<endl;
	for (it2= colLookUpTable.begin(); it2!=colLookUpTable.end(); it2++) {
		cout<<it2->first<<", "<<it2->second<<endl;
	}
#endif



	if (it != colLookUpTable.end()) {
#ifdef DBG
		cout<<"gene "<<it->first<<"'s index="<<it->second<<endl;
#endif
		return getColAsIntByOrder(it->second, order);
	} else {
		cerr<<"Error in table::getColAsIntByOrder, The table doesn't have column "<<name<<"."<<endl;
		vector<int> ret;
		return ret; //return an empty vector
	}
}


vector<int> table::getStatusByOrder(string st, vector<int> order) {
	map<string,int>::iterator it = colLookUpTable.find(st);
	if (it != colLookUpTable.end()) {
		return getStatusByOrder(it->second, order);
	} else {
		cerr<<"Error in table::getStatusByOrder, The table doesn't have column "<<st<<"."<<endl;
		vector<int> ret;
		return ret; //return an empty vector
	}
}

vector<int> table::getStatusByOrder(int j, vector<int> order) {
//  DECEASED  LIVING 
	int r = rows();
	vector<int> ret;
	for (int i=0; i<r; i++) {
//		cout<<"In table::getStatusByOrder, i="<<i<<endl;
		int status;
		string status_str = data[ order[i] ][j];
/*		cout<<"status_str = *"<<status_str<<"*"<<endl;
		cout<<status_str.compare("DECEASED")<<endl;
		cout<<status_str.compare("1")<<endl;
		cout<<status_str.compare("LIVING")<<endl;
		cout<<status_str.compare("0")<<endl;
*/
		if (status_str.compare("DECEASED")==0 || status_str.compare("1")==0)
			status = 1;
		else if (status_str.compare("LIVING")==0 || status_str.compare("0")==0)
			status = 0;
		//else
			//error
		ret.push_back( status );
	}
/*	cout << "ret = ";
	for (int i=0; i<ret.size(); i++){
		cout << ret[i] <<",";
	}
	cout << endl;
*/
	return ret;
}

int parseHeader(string line, vector<string> *out) {

#ifdef DBG
	cout<<"in parseHeader, line="<<line<<endl;
#endif

	int i=0, n=line.length();
	string tmp_str;
	tmp_str.clear();

	for (int i=0; i<n; i++) {

		//cout<<int(' ')<<" "<<int('\t')<<endl;
		//cout<<"line["<<i<<"]="<<(int)line[i]<<endl;
		if (line[i] != ' ' && line[i] != '\t') {
			//cout<<"line["<<i<<"]="<<line[i]<<endl;
			//cout<<"tmp_str="<<tmp_str<<endl<<endl;
			tmp_str+=line[i];

			if(i==n-1)
				out->push_back(tmp_str);

#ifdef DBG
			cout<<tmp_str<<endl;
			cout<<"------------------"<<endl;
#endif
		} else { 
			//cout<<"in else"<<endl;
			//cout<<"out->size()="<<out->size()<<endl;
			if (tmp_str.length() != 0) {
				out->push_back(tmp_str);
#ifdef DBG
				cout<<"pushed "<<tmp_str<<endl;
#endif
				tmp_str.clear();
			}
		}
		//cout<<endl<<endl;
	}

	n = out->size();
#ifdef DBG
	//cout<<"in parseHeader: n="<<n<<endl;
	for (int i=0; i<n; i++) {
		cout<<out->at(i)<<" ";
	}
	cout<<endl;
#endif

	return 0;
}

//use in gather sorted vector
bool compareS2B(struct order_pair a, struct order_pair b) {
	return a.entry<b.entry;
}
bool compareB2S(struct order_pair a, struct order_pair b) {
	return a.entry>b.entry;
}


void disp_vec_int(vector<int> v, string name) {
	int n = v.size();
	cout<<name<<" = "<<endl;
	cout<<"[";
	for (int i=0; i<n; i++) {
		cout<<v[i]<<" ";
	}
	cout<<"]"<<endl;
}


void disp_vec_string(vector<string> v, string name) {
	int n = v.size();
	cout<<name<<" = "<<endl;
	cout<<"[";
	for (int i=0; i<n; i++) {
		cout<<v[i]<<" ";
	}
	cout<<"]"<<endl;
}

