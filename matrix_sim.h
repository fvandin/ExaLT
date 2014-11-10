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

#ifndef MATRIX_SIM_H
#define MATRIX_SIM_H


#include<iostream>

#include<vector>

using namespace std;

class matrix_sim 
{
	private:
		vector < vector<double> > data;
		int row;
		int col;
	public:
		matrix_sim(int r, int c) ;
		void display(); 
		void joint_display( matrix_sim V,  double n, double n1, double e1);
		void zeros(int r, int c) ;
		void set_all(double d) ;
		void times(double d) ;
		vector<double> get_col(int c) ;
		matrix_sim get_col_as_row_matrix(int c) ;
		void set_element(int i, int j, double d) ;
		double get_element(int i, int j) ;
		int n_col();
		vector<double> get_row(int c);

};


#endif
