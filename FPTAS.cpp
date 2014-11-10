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
#include <cstdlib>


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "matrix_sim.h"
#include "FPTAS.h"

//#define DBG


using namespace std;

//%PRECISION on the values v_k; their approx. should not be too bad
double precision = 1e-10;

void display_vec(vector<double> v) {
	int size = v.size();
	cout<<"[ ";
	for (int i=0; i<size; i++) {
		cout<<v[i]<<" ";
	}
	cout<<"]"<<endl;
}


int vec_sum(vector<int> v) {
	int n = v.size();
	int sum = 0;
	for (int i=0; i<n; i++) {
		sum += v[i];
	}
	return sum;
}

int intervalBound(matrix_sim vtr, double value, double precision) {
	//%Searches the input vector to find l such that
	//%vector_{l-1}< value <=vector_(l)
	//%precision is used to define the accuracy of equality to the extremes of
	//%the vector; if this is not considered, problems may arise due to precision
	//%errors
	//
	int n_v = vtr.n_col();
	int l = 0;
	while ( ((l+1) <= n_v) && (double)value <= (double) vtr.get_element(1,l+1)){
		l++;
	}
	//adjust based on accuracy
	if ( ((l+1) <= n_v) && (fabs((double)value - vtr.get_element(1,l+1))<precision)){
		l++;
	}
	
	if (l==0){
		return l-1;
	}
	return l;
}



int intervalBound_secondSide(matrix_sim vtr, double value, double precision) {
	//%Searches the input vector to find l such that
	//%vector_{l-1}<=value<=vector_(l)
	//%precision is used to define the accuracy of equality to the extremes of
	//%the vector; if this is not considered, problems may arise due to precision
	//%errors
	//

	int n_v = vtr.n_col();
	int l = 0;
	while (  ((l+1) <= n_v) && (double)value >= (double) vtr.get_element(1,l+1)){
		l++;
	}
	//adjust based on accuracy
	if (((l+1) <= n_v) && (fabs((double)value - vtr.get_element(1,l+1))<precision)){
		l++;
	}

	if (l == 0){
		return l-1;
	}

	return l;

}



void array_max_min(double *a, int size, double *max_a, double *min_a) {
	*max_a = -DBL_MAX;
	*min_a =  DBL_MAX;
	for (int i=0; i<size; ++i) {
		if (a[i]>*max_a)
			*max_a = a[i];
		if (a[i]<*min_a)
			*min_a = a[i];
	}
}


void vector_max_min(vector<double> a, double *max_a, double *min_a) {
	*max_a = -DBL_MAX;
	*min_a =  DBL_MAX;
	int size = a.size();
	for (int i=0; i<size; ++i) {
		if (a[i]>*max_a)
			*max_a = a[i];
		if (a[i]<*min_a)
			*min_a = a[i];
	}
}

void display_matrix(vector< vector<double> > v) {
	int row = v.size();
	if (row == 0)
		return;
	int col = v[0].size();
	if (col==0)
		return;
	int i,j;
	for (i=0;i<row;i++) {
		for (j=0;j<col;j++) {
			cout<<v[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;
}


void display_matrix_array(double **v, int row, int col) {
	int i,j;
	for (i=0;i<row;i++) {
		for (j=0;j<col;j++) {
			cout<<v[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;
}


void display_array(double *a, int size) {
	int i;
	cout<<"[ ";
	for (i=0;i<size;i++) {
		cout<<a[i]<<"  ";
	}
	cout<<"]"<<endl;;
}

int FPTAS(vector<int> x_v, vector<int> c_v, double e, double &p_value, matrix_sim &v, matrix_sim &p) {

	//check to make sure x_v and c_v has the same size
	if ( x_v.size() != c_v.size() ) {
		cerr<<"Error in FPTAS: input x and c must have the same size, x.size()="<<x_v.size()<<" while c.size()="<<c_v.size()<<"."<<endl;
		fflush(stderr);
		return -1;
	}

	/*
	   %
	   % A Fully Polynomial Time Approximation Scheme for the Unconditional LogRank Test
	   %
	   % INPUT
	   % x: row vector of n events, n1 of them happening in G1. x_j=1 if j-th event
	   % comes from G1, 0 if from G0.
	   % c: row vector of m censored events. c_j=0 if j-th event is censored, 0
	   % otherwise
	   % e: approximation, 0<e<=1
	   %
	   % OUTPUT
	   %
	   */

	//TODO: change the approximation of (n \choose n_1) used to compute the minimum probability (if n_1 >=3); that will reduce the maximum k used

	//% %Initialization
	int n  = x_v.size(); //%number of events
	int n1 = vec_sum(x_v);   //%number of G0 events
	//FABIO 08-04: new way to compute e1 BEGIN
//	double e1=e/(2*n);  //%approximation factor
	double e1 = 1.0 - pow( 1.0+e, -1.0/((double)n) );
	//change x_v, c_v to x, c
	matrix_sim x(1, n);
	matrix_sim c(1, n);
	for (int i=0; i<n; i++) {
		x.set_element(1, i+1, x_v[i]);
		c.set_element(1, i+1, c_v[i]);
	}
	double tmp_dbl = (double)n1*log( (double)n )/log(1+e1);
	int k_max = tmp_dbl - (int)tmp_dbl > 0? (int)tmp_dbl+1 :  (int)tmp_dbl; //%maximum number of k
	int k_max_index = k_max + 1;         //%+1 since Matlab starts in position 1
	int n1_index = n1 + 1;
	//cout << "n1_index: " << n1_index << endl;
	v.zeros(k_max_index, n1_index);
	v.set_all(DBL_MAX);
	p.zeros(k_max_index, n1_index);
	p.set_all(0.0);
	//%values of probabilities to check for k
	matrix_sim p_k(1,1);
	p_k.zeros(1, k_max+1);
	p_k.set_all(0.0);
	for (int k=0; k<=k_max; k++) 
		p_k.set_element(1, k+1,   ( pow(n,-n1)  ) * (  pow((1 + e1),k) ) );
	v.set_element(k_max_index,1,0);  // %v(k_max_index,1)=0;
	p.set_element(k_max_index,1,1);
	matrix_sim p_prime(1,1);
	matrix_sim v_prime(1,1);
	int r_index, k_index;
	double u_prime_prime, u_prime, v1, v2;
	double p1, p2;
	int l;
	
	//%t=current position in the string
	for (int t=1; t<=n; t++) {
		p_prime.zeros(2*k_max_index,n1_index);
		p_prime.set_all(0.0);
		v_prime.zeros(2*k_max_index,n1_index);
		v_prime.set_all(DBL_MAX);
		//%iterate over all values k to identify the values of v(k,r) that are
		//%\neq -Inf
		int min_n1_t = min(n1,t);
		for (int r=0; r<=min_n1_t ; r++) {
			r_index = r + 1;
			for (int k=0; k<=k_max; k++) {
				k_index = k + 1;
				if (v.get_element(k_index, r_index) != DBL_MAX) {
					//check if same element already computed!
					if ((k_index == 1) || v.get_element(k_index, r_index) != v.get_element(k_index-1, r_index)){
						if (c.get_element(1,t) == 1) {
							if (r < n1) {
								/*
								%new value of the statistic if a one is found in
								%position t; this can happen only if not all 1's
								%have been placed, i.e. r<n1
								*/
								
								v1=v.get_element(k_index,r_index)+1-((double)n1-r)/((double)n-t+1);
								/*
								%value v_1 corresponds to have r+1 events; the
								%value can be obtained: in two ways
								%1) from having before r+1
								%events and having now chosen a 0 (that happens 
								%with probability (1 - (n1 -(r+1))/n_t))
								%2) from having r events before and having now
								%chosen a 1, that happens with probability
								%(n1-r)/n_t
								%case 1): before r+1 events, and 0 chosen
								*/
								u_prime=v1 + ((double)n1-(r+1))/((double)n-t+1);
								l = intervalBound(v.get_col_as_row_matrix(r_index+1), u_prime, precision); //%14
								//if l~=-1
								if (l != -1){
									p1=p.get_element(l,r_index+1)*(1 - ((double)n1-(r+1))/((double)n-t+1));
								} else
									p1 = 0;
								//%case 2): before r events, and 1 chosen
								u_prime_prime = v1-(1-(double)(n1-r)/(n-t+1));
								l = intervalBound(v.get_col_as_row_matrix(r_index), u_prime_prime, precision);
								if ( l!=-1 ) {
									p2=p.get_element(l,r_index)*((double)n1-r)/((double)n-t+1);
								}
								else
									p2 = 0;
								p_prime.set_element((2*k_index-1),(r_index+1), p1+p2 );
								v_prime.set_element((2*k_index-1),(r_index+1), v1 );                  
							} 
							//%new value of the statistic if a 0 is found in position
							//%t
							v2=v.get_element(k_index, r_index) - ((double)n1-r)/((double)n-t+1); //%event in G0
							/*
							%value v_2 corresponds to have r events and therefore a 0 in current position; the
							%value can be obtained: in two ways
							%1) from having before r
							%events and having now chosen a 0 (that happens 
							%with probability (1 - (n1 -r)/n_t))
							%2) from having r-1 events before and having now
							%chosen a 1, that happens with probability
							%(n1-(r-1))/n_t
							%case 1)
							*/
							u_prime=v2+((double)n1-r)/((double)n-t+1);
							l = intervalBound(v.get_col_as_row_matrix(r_index), u_prime, precision);
							if (l!=-1) {
								p1=p.get_element(l,r_index)*(1-((double)n1-r)/((double)n-t+1));
							} else {
								p1=0;
							}
							//%case 2; possible only if r is > 0
							if (r>0) {
								u_prime_prime = v2-(1-((double)n1-(r-1))/((double)n-t+1));
								l=intervalBound(v.get_col_as_row_matrix(r_index-1), u_prime_prime, precision);
								if (l != -1) {
									p2=p.get_element(l,r_index-1)*((double)n1-(r-1))/((double)n-t+1);
								} else
									p2=0;
							} else
								p2=0;
							p_prime.set_element((2*k_index),r_index,p1+p2);
							v_prime.set_element(2*k_index,r_index,v2);
						} else  {  //%there is censoring in position t (c(t) == 0)
							//%case 1 found in the current position
							if ( r<n1 ) {
								u_prime=v.get_element(k_index,r_index);
								l = intervalBound(v.get_col_as_row_matrix(r_index+1), u_prime, precision);
								if (l != -1) {
									p1=p.get_element(l,r_index+1)*(1-((double)n1-(r+1))/((double)n-t+1));
								} else {
									p1 = 0;
								}
								l=intervalBound(v.get_col_as_row_matrix(r_index), u_prime, precision);
								if (l!=-1) {
									p2=p.get_element(l,r_index)*((double)n1-r)/((double)n-t+1);
								} else
									p2 = 0;
								p_prime.set_element((2*k_index-1),r_index+1, p1+p2);
								v_prime.set_element((2*k_index-1),r_index+1, u_prime);
							} 
							//%case 0 found in the current position
							u_prime=v.get_element(k_index,r_index);
							l=intervalBound(v.get_col_as_row_matrix(r_index),u_prime, precision);
							if (l != -1) {
								p1=p.get_element(l,r_index)*(1-((double)n1-r)/((double)n-t+1));
							} else
								p1 = 0;

							if (r>0) {
								l=intervalBound(v.get_col_as_row_matrix(r_index-1), u_prime, precision);
								if (l != -1) {
									p2=p.get_element(l,r_index-1)*((double)n1-(r-1))/((double)n-t+1);
								} else
									p2=0;
							} else
								p2=0;
							p_prime.set_element((2*k_index),r_index, p1+p2);
							v_prime.set_element((2*k_index),r_index, u_prime);
						}
					}
				} 
			} //k loop 
		} //r loop
		//%build the new tables
		v.zeros(k_max_index,n1_index); //%logrank statistics
		v.set_all(DBL_MAX);

		p.zeros(k_max_index, n1_index);
		p.set_all(-1.0);
		for (int r=0; r<=n1; r++ ) { 
			r_index = r + 1;
			//for k=0:1:k_max
			for (int k=0;k<=k_max;k++) {
				k_index = k + 1;
				double max_v = -DBL_MAX;
				double max_prob = -DBL_MAX;
//				int index_max_v=0;
				double lower_b = pow(n,-n1) *pow(1+e1,k-1);
				if (k == 0){
					lower_b = 0.0;
				}
				int k_max_index_2times = 2*k_max_index;
				for (int i=1; i<=k_max_index_2times; i++) {
					if ( p_prime.get_element(i,r_index)<=pow(n,-n1) * pow(1+e1,k) &&  p_prime.get_element(i,r_index)>lower_b)
					{
						if (v_prime.get_element(i,r_index)>max_v){
							max_v=v_prime.get_element(i,r_index);
						}
						if ( p_prime.get_element(i,r_index) > max_prob){
							max_prob = p_prime.get_element(i,r_index);
						}
//						index_max_v=i;
					}
				}
				if (max_v!=-DBL_MAX) {
					v.set_element(k_index,r_index, max_v);
					//%min_v
					p.set_element(k_index,r_index, max_prob);
				} 
			} 
			//now repeat elements whenever needed
			double curr_v = DBL_MAX;
			double curr_prob = 0.0;
			for (int k=0;k<=k_max;k++){
				k_index = k+1;
				if (p.get_element(k_index,r_index) == -1.0){
					p.set_element(k_index,r_index, curr_prob);
					v.set_element(k_index,r_index, curr_v);
				}
				else{
					curr_prob = p.get_element(k_index,r_index);
					curr_v = v.get_element(k_index, r_index);
				}
			}
		}
	}//t loop
	//% P-value calculation
	double V = 0;
	for (int j=1;j<=n;j++) {
		double sum_x_0_j_1 = 0;
		for (int ii=1; ii<=j-1; ii++) {
			sum_x_0_j_1 += x.get_element(1, ii);
		}
		V = V + c.get_element(1,j) * (  (double)x.get_element(1,j) -(n1-sum_x_0_j_1 )/((double)n-j+1));
	}

	V=fabs(V);
	double max_v_column_n1_index = -DBL_MAX;
	double min_v_column_n1_index = DBL_MAX;
	vector_max_min(v.get_col(n1_index), &max_v_column_n1_index, &min_v_column_n1_index);

	//void array_max_min(vector<double> a, double *max_a, double *min_a) 
	//double p[k_max_index][n1_index];
	//double p_value;
	double min_p_column_n1_index, max_p_column_n1_index;
	vector_max_min(p.get_col(n1_index), &max_p_column_n1_index, &min_p_column_n1_index);
	if ( V>max_v_column_n1_index ) {
		p_value = 0;
	}
	else if ( V<=min_v_column_n1_index ) {
		p_value = max_p_column_n1_index;
	} else {
		int l = intervalBound(v.get_col_as_row_matrix(n1_index), V, precision);
		p_value=p.get_element(l,n1_index);
	}

	return 0;
}




void FPTAS_secondSide(vector<int> x_v, vector<int> c_v, double e, double &p_value, matrix_sim &v, matrix_sim &p) {

	//%PRECISION on the values v_k; their approx. should not be too bad
	double precision = 1e-10;
	/*
	   %
	   % A Fully Polynomial Time Approximation Scheme for the Unconditional LogRank Test
	   %
	   % INPUT
	   % x: row vector of n events, n1 of them happening in G1. x_j=1 if j-th event
	   % comes from G1, 0 if from G0.
	   % c: row vector of m censored events. c_j=0 if j-th event is censored, 0
	   % otherwise
	   % e: approximation, 0<e<=1
	   %
	   % OUTPUT
	   %
	   */


	//% %Initialization
	int n  = x_v.size(); //%number of events
	int n1 = vec_sum(x_v);   //%number of G0 events
        double e1 = 1.0 - pow( 1.0+e, -1.0/((double)n) ); //approximation factor

	//change x_v, c_v to x, c
	matrix_sim x(1, n);
	matrix_sim c(1, n);
	for (int i=0; i<n; i++) {
		x.set_element(1, i+1, x_v[i]);
		c.set_element(1, i+1, c_v[i]);
	}



	double tmp_dbl = (double)n1*log( (double)n )/log(1+e1);
	int k_max = tmp_dbl - (int)tmp_dbl > 0? (int)tmp_dbl+1 :  (int)tmp_dbl; //%maximum number of k
	int k_max_index = k_max + 1;         //%+1 since Matlab starts in position 1
	int n1_index = n1 + 1;

	v.zeros(k_max_index, n1_index);
	v.set_all(-DBL_MAX);
	p.zeros(k_max_index, n1_index);
	p.set_all(0.0);


	//%values of probabilities to check for k
	matrix_sim p_k(1,1);
	p_k.zeros(1, k_max+1);
	p_k.set_all(0.0);
	for (int k=0; k<=k_max; k++) 
		p_k.set_element(1, k+1,   ( pow(n,-n1)  ) * (  pow((1 + e1),k) ) );

	v.set_element(k_max_index,1,0);  // %v(k_max_index,1)=0;
	p.set_element(k_max_index,1,1);

	matrix_sim p_prime(1,1);
	matrix_sim v_prime(1,1);
	int r_index, k_index;
	double u_prime_prime, u_prime, v1, v2;
	double p1, p2;
	int l;

	//%t=current position in the string
	for (int t=1; t<=n; t++) {

		p_prime.zeros(2*k_max_index,n1_index);
		p_prime.set_all(0.0);
		v_prime.zeros(2*k_max_index,n1_index);
		v_prime.set_all(-DBL_MAX);

		//%iterate over all values k to identify the values of v(k,r) that are
		//%\neq -Inf
		int min_n1_t = min(n1,t);

		for (int r=0; r<=min_n1_t ; r++) {
			r_index = r + 1;
			for (int k=0; k<=k_max; k++) {
				k_index = k + 1;
				if (v.get_element(k_index, r_index) != -DBL_MAX) {
					if (c.get_element(1,t) == 1) {
						if (r < n1) {
							/*
							   %new value of the statistic if a one is found in
							   %position t; this can happen only if not all 1's
							   %have been placed, i.e. r<n1
							   */
							v1=v.get_element(k_index,r_index)+1-((double)n1-r)/((double)n-t+1);
							/*

							   %value v_1 corresponds to have r+1 events; the
							   %value can be obtained: in two ways
							   %1) from having before r+1
							   %events and having now chosen a 0 (that happens 
							   %with probability (1 - (n1 -(r+1))/n_t))
							   %2) from having r events before and having now
							   %chosen a 1, that happens with probability
							   %(n1-r)/n_t
							   %case 1): before r+1 events, and 0 chosen
							   */
							u_prime=v1 + ((double)n1-(r+1))/((double)n-t+1);
							l = intervalBound_secondSide(v.get_col_as_row_matrix(r_index+1), u_prime, precision); //%14
							//if l~=-1
							if (l != -1){

								p1=p.get_element(l,r_index+1)*(1 - ((double)n1-(r+1))/((double)n-t+1));
							} else
								p1 = 0;
							//%case 2): before r events, and 1 chosen
							u_prime_prime = v1-(1-(double)(n1-r)/(n-t+1));
							l = intervalBound_secondSide(v.get_col_as_row_matrix(r_index), u_prime_prime, precision);
							if ( l!=-1) {
								p2=p.get_element(l,r_index)*((double)n1-r)/((double)n-t+1);
							}
							else
								p2 = 0;
							p_prime.set_element((2*k_index-1),(r_index+1), p1+p2 );
							v_prime.set_element((2*k_index-1),(r_index+1), v1 );                  
						} 
						//%new value of the statistic if a 0 is found in position
						//%t
						v2=v.get_element(k_index, r_index) - ((double)n1-r)/((double)n-t+1); //%event in G0
						/*
						   %value v_2 corresponds to have r events and therefore a 0 in current position; the
						   %value can be obtained: in two ways
						   %1) from having before r
						   %events and having now chosen a 0 (that happens 
						   %with probability (1 - (n1 -r)/n_t))
						   %2) from having r-1 events before and having now
						   %chosen a 1, that happens with probability
						   %(n1-(r-1))/n_t
						   %case 1)
						   */
						u_prime=v2+((double)n1-r)/((double)n-t+1);
						l = intervalBound_secondSide(v.get_col_as_row_matrix(r_index), u_prime, precision);
						if (l!=-1) {
							p1=p.get_element(l,r_index)*(1-((double)n1-r)/((double)n-t+1));
						} else {
							p1=0;
						}
						//%case 2; possible only if r is > 0
						if (r>0) {
							u_prime_prime = v2-(1-((double)n1-(r-1))/((double)n-t+1));
							l=intervalBound_secondSide(v.get_col_as_row_matrix(r_index-1), u_prime_prime, precision);
							if (l != -1) {
								p2=p.get_element(l,r_index-1)*((double)n1-(r-1))/((double)n-t+1);
							} else
								p2=0;
						} else
							p2=0;
						p_prime.set_element((2*k_index),r_index,p1+p2);
						v_prime.set_element(2*k_index,r_index,v2);
					} else  {  //%there is censoring in position t (c(t) == 0)
						//%case 1 found in the current position
						if ( r<n1 ) {
							u_prime=v.get_element(k_index,r_index);
							l = intervalBound_secondSide(v.get_col_as_row_matrix(r_index+1), u_prime, precision);
							if (l != -1) {
								p1=p.get_element(l,r_index+1)*(1-((double)n1-(r+1))/((double)n-t+1));
							} else {
								p1 = 0;
							}
							l=intervalBound_secondSide(v.get_col_as_row_matrix(r_index), u_prime, precision);
							if (l!=-1) {
								p2=p.get_element(l,r_index)*((double)n1-r)/((double)n-t+1);
							} else
								p2 = 0;
							p_prime.set_element((2*k_index-1),r_index+1, p1+p2);
							v_prime.set_element((2*k_index-1),r_index+1, u_prime);
						} 
						u_prime=v.get_element(k_index,r_index);
						l=intervalBound_secondSide(v.get_col_as_row_matrix(r_index),u_prime, precision);
						if (l != -1) {
							p1=p.get_element(l,r_index)*(1-((double)n1-r)/((double)n-t+1));
						} else
							p1 = 0;
						if (r>0) {
							l=intervalBound_secondSide(v.get_col_as_row_matrix(r_index-1), u_prime, precision);
							if (l != -1) {
								p2=p.get_element(l,r_index-1)*((double)n1-(r-1))/((double)n-t+1);
							} else
								p2=0;
						} else
							p2=0;
						p_prime.set_element((2*k_index),r_index, p1+p2);
						v_prime.set_element((2*k_index),r_index, u_prime);
					} 
				} 
			} //k loop 
		} //r loop

		//%build the new tables

		v.zeros(k_max_index,n1_index); //%logrank statistics
		v.set_all(-DBL_MAX);

		p.zeros(k_max_index, n1_index);
		p.set_all(-1.0);
		
		//for r=0:1:min(n1,t)    %OR n1 as in the pdf?
		for (int r=0; r<=n1; r++ ) { //   %OR n1 as in the pdf?
			r_index = r + 1;
			//for k=0:1:k_max
			for (int k=0;k<=k_max;k++) {
				
				double lower_b = pow(n,-n1) *pow(1+e1,k-1);
				if (k == 0){
					lower_b = 0.0;
				}
				
				k_index = k + 1;

				double min_v = DBL_MAX;
//				int index_max_v=0;
				double max_prob = -DBL_MAX;
				int k_max_index_2times = 2*k_max_index;
				for (int i=1; i<=k_max_index_2times; i++) {
					if ( (p_prime.get_element(i,r_index)<=(pow(n,-n1) * pow(1+e1,k)) ) && p_prime.get_element(i,r_index) > lower_b )
					{
						if (v_prime.get_element(i,r_index)<min_v){
							min_v=v_prime.get_element(i,r_index);
						}
						if ( p_prime.get_element(i,r_index) > max_prob){
							max_prob = p_prime.get_element(i,r_index);
						}
					}
				} 
				if (min_v!=DBL_MAX) {
					v.set_element(k_index,r_index, min_v);
					//%min_v
					p.set_element(k_index,r_index, max_prob);
				} 
			} 
			//now repeat elements whenever needed
			double curr_v = -DBL_MAX;
			double curr_prob = 0.0;
			for (int k=0;k<=k_max;k++){
				k_index = k+1;
				if (p.get_element(k_index,r_index) == -1.0){
					p.set_element(k_index,r_index, curr_prob);
					v.set_element(k_index,r_index, curr_v);
				}
				else{
					curr_prob = p.get_element(k_index,r_index);
					curr_v = v.get_element(k_index, r_index);
				}
			}

		}
	}//t loop

	//% P-value calculation
	double V = 0;
	double seen = 0;
	for (int j=1;j<=n;j++) {

		V = V + c.get_element(1,j) * (  (double)x.get_element(1,j) -(n1-seen)/((double)n-j+1));
		seen += (double)x.get_element(1,j);
	}

	V=-fabs(V);

	double max_v_column_n1_index = -DBL_MAX;
	double min_v_column_n1_index = DBL_MAX;
	vector_max_min(v.get_col(n1_index), &max_v_column_n1_index, &min_v_column_n1_index);

	//void array_max_min(vector<double> a, double *max_a, double *min_a) 
	//double p[k_max_index][n1_index];
	//double p_value;
	double min_p_column_n1_index, max_p_column_n1_index;
	vector_max_min(p.get_col(n1_index), &max_p_column_n1_index, &min_p_column_n1_index);

	if ( V>max_v_column_n1_index ) {
		//find the minimum p greater than 0
		vector<double> tmp = p.get_col(n1_index);
		int size = tmp.size();
		double min_p = DBL_MAX;
		for (int i=0; i<size; ++i) {
			if ((tmp[i]<min_p) && (tmp[i]!=0)){
				min_p = tmp[i];
			}
		}
		p_value = min_p;
	}
	else if ( V<=min_v_column_n1_index ) {
		p_value=p.get_element(k_max_index-1,n1_index);
	} else {

		int l = intervalBound_secondSide(v.get_col_as_row_matrix(n1_index), V, precision);
		//l
		p_value=p.get_element(l,n1_index);
	}
}



