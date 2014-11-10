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


#ifndef FPTAS_H
#define FPTAS_H

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

//#define DBG


using namespace std;

void display_vec(vector<double> v);


int vec_sum(vector<int> v); 

int intervalBound(matrix_sim vtr, double value, double precision); 


void array_max_min(double *a, int size, double *max_a, double *min_a); 


void vector_max_min(vector<double> a, double *max_a, double *min_a); 

void display_matrix(vector< vector<double> > v); 


void display_matrix_array(double **v, int row, int col); 


void display_array(double *a, int size); 


int FPTAS(vector<int> x_v, vector<int> c_v, double e, double &p_value, matrix_sim &v, matrix_sim &p); 


void FPTAS_secondSide(vector<int> x_v, vector<int> c_v, double e, double &p_value, matrix_sim &v, matrix_sim &p); 





#endif
