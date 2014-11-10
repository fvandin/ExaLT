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


#include "matrix_sim.h"
#include <stdio.h>
#include <math.h>

matrix_sim::matrix_sim(int r, int c) 
{
	row = r;
	col = c;
	data.clear();
	int i,j;
	data.reserve(row);
	for (i=0; i<row; i++) 
	{
		vector<double> tmp;
		for (j=0; j<col; j++) 
		{
			tmp.reserve(col);
			tmp.push_back(0.0);
		}
		data.push_back(tmp);
	}

}


void matrix_sim::display() 
{
	int i,j;
	for (i=0; i<row; i++) 
	{
		for (j=0; j<col; j++) 
		{
			cout<<data[i][j]<<"  ";
		}
		cout<<endl;
	}
}

void matrix_sim::joint_display( matrix_sim V, double n, double n1, double e1 ) 
{
	int i,j;
	for (i=0; i<row; i++) 
	{
		cout <<"[" << ( pow(n,-n1)  ) * (  pow((1 + e1),i) )<<"]\t";
		for (j=0; j<col; j++) 
		{
			cout<<"("<<data[i][j]<<", "<< V.get_element(i+1,j+1)<<")\t";
		}
		cout<<endl;
	}
}

void matrix_sim::zeros(int r, int c) 
{
	row = r;
	col = c;
	data.clear();
	int i,j;
	data.reserve(r);
	for (i=0; i<row; i++) 
	{
		vector<double> tmp;
		tmp.reserve(col);
		for (j=0; j<col; j++) 
		{
			tmp.push_back(0.0);
		}
		data.push_back(tmp);
	}

}

void matrix_sim::set_all(double d) 
{
	int i,j;
	for (i=0; i<row; i++) 
	{
		for (j=0; j<col; j++) 
		{
			data[i][j] = d;
		}
	}
}

void matrix_sim::times(double d) 
{
	int i,j;
	for (i=0; i<row; i++) 
	{
		for (j=0; j<col; j++) 
		{
			data[i][j] *= d;
		}
	}
}

vector<double> matrix_sim::get_col(int c) 
{
	vector<double> ret;
	int i,j;
	for (i=0; i<row; i++) 
	{
			ret.push_back(data[i][c-1]);
	}
	return ret;
}


matrix_sim matrix_sim::get_col_as_row_matrix(int c) 
{
	matrix_sim ret(1, row);
	int i,j;
	//cout<<"row="<<endl;
	for (i=1; i<=row; i++) 
	{
		//cout<<"i="<<i<<endl;
		ret.set_element(1, i,  (data[i-1][c-1]) );
		//cout<<"i="<<i<<" good"<<endl;
	}
	return ret;
}

void matrix_sim::set_element(int i, int j, double d) 
{
	if (i<1) 
	{
		cerr<<"Error, i<1"<<endl;
		return;
	}
	if (j<1) 
	{
		cerr<<"Error, j<1"<<endl;
		return;
	}



	/*
	cout<<data.size
	cout<<"debug: i="<<i<<", j="<<j<<endl; */
	data[i-1][j-1] = d;
}

double matrix_sim::get_element(int i, int j) 
{
	return data[i-1][j-1];
}


int matrix_sim::n_col()
{
	return col;
}

vector<double> matrix_sim::get_row(int c) 
{
	vector<double> ret;
	int i;
	for (i=0; i<col; i++) 
	{
			ret.push_back(data[c-1][i]);
	}
	return ret;
}