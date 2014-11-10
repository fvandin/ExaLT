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

#include<Rcpp.h>


#include<table.h>
#include<FPTAS.h>
#include<logrank_exact.h>
#include<MC_logrank.h>

#include <string>
#include <pthread.h>


#define min(a,b) a<=b?a:b

using namespace std;


typedef struct {
	vector<int> *x ;
	vector<int> *c ;
	double *e;
	double *p_value;
	matrix_sim *v;
	matrix_sim *p;
} pt_data;



void *FPTAS1(void* t) {
	pt_data *data1 = (pt_data *)t;
	FPTAS( *(data1->x), *(data1->c), *(data1->e), *(data1->p_value), *(data1->v), *(data1->p) );
}

void *FPTAS2(void* t) {
	pt_data *data1 = (pt_data *)t;
	FPTAS_secondSide( *(data1->x), *(data1->c), *(data1->e), *(data1->p_value), *(data1->v), *(data1->p) );
}


//RcppExport SEXP R_FPTAS(SEXP table_file, SEXP gene) {
RcppExport SEXP R_ExaLT(SEXP table_file, SEXP gene) {


	//Rcpp::StringVector table_name(table_file);
	string name = Rcpp::as<std::string>(table_file);
	string gene_str = Rcpp::as<std::string>(gene);


	//build table
	table t((char*)name.c_str());

	//get x,c by sorted time order 
	vector<int> order = t.sort_order(2, SMALL_TO_LARGE);  //time col index is 2
	vector<int> c = t.getStatusByOrder("STATUS", order);
	vector<int> c2 = c;
	vector<int> x = t.getColAsIntByOrder((char*)gene_str.c_str(), order);
	vector<int> x2 = x;

	//call the function
	double p_value;
	double p_value2;

	//run exhaustive if the number of combinations is not too large
	unsigned long long n_ones = 0;
	for (int i=0;i<x.size(); i++){
		if (x[i]==1) {
			n_ones += 1;
		}
	}
	cout <<"number of samples in group 0: "<<x.size()-n_ones<<endl;
	cout <<"number of samples in group 1: "<<n_ones<<endl;
	unsigned long long ncomb = bcoeff(x.size() , n_ones);
	if (ncomb <= 100000000){
		cout<<"Computing exact p-value..."<<endl;
		double p_value_exh;
		double p_value_exh_left;
		double p_value_exh_right;
		exhaustive( x, c, p_value_exh, p_value_exh_left, p_value_exh_right);
		cout <<"p-value (left) = " << p_value_exh_left << endl;
		cout <<"p-value (right) = " << p_value_exh_right << endl;
		cout<<"p-value = "<<p_value_exh<<endl;
		Rcpp::NumericVector pvalue_exh(1);
                pvalue_exh[0] = p_value_exh;
                return(pvalue_exh);
}

	//run MC to estimate p-value first
	double p_value_MC;
	double p_value_MC_left;
	double p_value_MC_right;
	unsigned long long number_iter = 10000000;
	cout<<"Running quick estimate of the p-value..."<<endl;
	MC_estimate( x, c, p_value_MC, p_value_MC_left, p_value_MC_right, number_iter );
	if (p_value_MC == 0.0 ){
		double min_pvalue_MC = 1.0/double(number_iter);
		p_value_MC_left,
		cout <<"the estimated p-value is < "<<min_pvalue_MC<<endl;
	}
	else{
		cout <<"[The left p-value is estimated to be close to: "<<p_value_MC_left << "]" <<endl;
		cout <<"[The right p-value is estimated to be close to: "<<p_value_MC_right << "]" <<endl;
		cout<<"the p-value is estimated to be close to: "<<p_value_MC<<endl;
	}

	//now check if more accurate approximation is needed

	char input;
	do{
		cout << "Do you want a controlled approximation? [Y/N]" << endl;
		cin >> input;
	} while ( (input != 'Y') && (input != 'N') );
	if (input == 'N') {
		Rcpp::NumericVector pvalue_MC(1);
                pvalue_MC[0] = p_value_MC;
                return(pvalue_MC);
	}
	else{

		cout << "Provide approximation factor (see README)"<<endl;
		double e;
		cin >> e;
		cout << "Obtaining controlled approximation..." << endl;

		matrix_sim v(1,1),p(1,1);
		matrix_sim v2(1,1),p2(1,1);


		pt_data data1;
		data1.x = &x;
		data1.c = &c;
		data1.e = &e;
		data1.p_value = &p_value;
		data1.v = &v;
		data1.p = &p;


		pt_data data2;
		data2.x = &x2;
		data2.c = &c2;
		data2.e = &e;
		data2.p_value = &p_value2;
		data2.v = &v2;
		data2.p = &p2;
		
		int rc;

		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE); 


		pthread_t thread_FPTAS1, thread_FPTAS2;
		rc = pthread_create(&thread_FPTAS1, &attr, FPTAS1, (void*)&data1);
		if (rc) {
			fprintf(stderr, "ERROR: return code from pthread_create() is %d\n", rc);
			exit(-1);
		}

		rc = pthread_create(&thread_FPTAS2, &attr, FPTAS2, (void*)&data2);
		if (rc) {
			fprintf(stderr, "ERROR: return code from pthread_create() is %d\n", rc);
			exit(-1);
		}

		//join the threads
		void* status;
		for (int i=0; i<2; i++) {
			if (i==0) {
				rc = pthread_join( thread_FPTAS1, &status);
			} else if (i==1) {
				rc = pthread_join( thread_FPTAS2, &status);
			}
			if (rc) {
				fprintf(stderr,  "return code from pthread_join() is %d\n", rc);
				exit(-1);
			}
		}
	
		cout<<"p-value (left) = "<<p_value2<<endl;
		cout<<"p-value (right)  = "<<p_value<<endl;

		p_value = min(p_value+p_value2, 1);
		cout<<"p-value = "<<p_value<<endl;
	

		Rcpp::NumericVector pvalue(1);
		pvalue[0] = p_value;

		return(pvalue);
	}
}


