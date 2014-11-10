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
#include "FPTAS.cpp"
#include "logrank_exact.h"
#include "MC_logrank.h"

#include <string.h>
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

//survival times vector
vector<float> survival_times;
//number of iterations for MC estimate
unsigned long long number_iter = 100000;

void *FPTAS1(void* t) {
	pt_data *data1 = (pt_data *)t;

	FPTAS( *(data1->x), *(data1->c), *(data1->e), *(data1->p_value), *(data1->v), *(data1->p) );
}

void *FPTAS2(void* t) {
	pt_data *data1 = (pt_data *)t;
	FPTAS_secondSide( *(data1->x), *(data1->c), *(data1->e), *(data1->p_value), *(data1->v), *(data1->p) );
}

//function used to derive permutation sorting for survival times
bool compare_surv_times(float a, float b) {
  return (survival_times[a] < survival_times[b]);
}

//RcppExport SEXP R_FPTAS(SEXP table_file, SEXP gene) {
int main (int argc, const char* argv[]) {

  	//TODO: check if table (with flag -table) is passed in input. If so, run automatically on all genes in the table.
	/*FORMAT TABLE:
	first line: header - first toke is patientID, second token is survival time, third is censoring info, then geneIDs, everything \t separated
	other lines: one entry per patient; patientID is string; survival time is float; censoring is 1 if death, is 0 otherwise; for geneID colums, 
		     is 0 if gene not mutated in the patient and
		     it is 1 otherwise.

	OUTPUT FORMAT:
	list of geneIDs and p-values and method used to compute p-value
	after -table, pass -p_threshold with threshold p-value for MC method -> if p-value estimated with
	MC > p_threshold then use that, otherwise use FPTAS; if not passed by user, use 0.05/num.genes
	- requires to estimate number of genes from table first
	*/
	string tmp = string( "-table" );
	if ( argc > 2 and tmp.compare( argv[2] ) == 0 ) {
	  cout << "RUNNING ON ENTIRE TABLE PASSED IN INPUT" << endl;
	  //order input params: table_file -table min_freq
	  //read input file
	  ifstream intable;
	  intable.open(argv[1],ifstream::in);
	  string line;
	  int line_c = 0;
	  int num_genes = 0;
	  int num_patients = 0;
	  //survival times vector
	  //vector <float> survival_times;
	  //censoring vector
	  vector <int> censoring_tmp;
	  //geneIDs (same order as in file)
	  vector <string> geneID;
	  if (intable.is_open()) {
	    while ( getline (intable,line) ){
	      //cout << line << endl;
	      char * line_char = strdup( line.c_str() );
	      int tok_counter = 0;
	      char *p = strtok(line_char, "\t");
	      while (p) {
		//printf ("Token: %s\n", p);
		if (line_c == 0 and tok_counter >=3){
		  num_genes++;
		  geneID.push_back(string(p));
		}
		if (line_c > 0 and tok_counter == 1){
		  survival_times.push_back(atof(p));
		}
		if (line_c > 0 and tok_counter == 2){
		  char* pEnd;
		  censoring_tmp.push_back(atoi(p));
		}
		tok_counter++;
		p = strtok(NULL, "\t");
	      }
	      if (line_c > 0){
		num_patients++;
	      }
	      line_c++;
	    }
	    intable.close();
	  }
	  cout << "Number of genes: " << num_genes << endl;
	  cout << "Number of patients: " << num_patients << endl;
 	  //read minimum frequency
	  int min_freq = atoi( argv[3] );
	  cout << "Minimum frequency: " << min_freq << endl;
	  /*
	  //print content survival_times
	  for (vector<float>::iterator it = survival_times.begin(); it != survival_times.end(); ++it)
	    cout << ' ' << *it;
	  cout << endl;
	  //print content censoring
	  for (vector<int>::iterator it = censoring_tmp.begin(); it != censoring_tmp.end(); ++it)
	    cout << ' ' << *it;
	  cout << endl;
	  */
	  //gene table
	  vector< vector<int> > gene_table( num_genes, vector<int>(0) );
	  line_c = 0;
	  intable.open(argv[1],ifstream::in);
	  if (intable.is_open()) {
	    while ( getline (intable,line) ){
	      char * line_char = strdup( line.c_str() );
	      int tok_counter = 0;
	      char *p = strtok(line_char, "\t");
	      int gene_index = 0;
	      while (p) {
		if (line_c > 0 and tok_counter >=3){
		  gene_table[gene_index].push_back(atoi(p));
		  gene_index++;
		}
		tok_counter++;
		p = strtok(NULL, "\t");
	      }
	      line_c++;
	    }
	    intable.close();
	  }
	  //print content gene_table
	  /*
	  for (int i=0; i<num_genes; i++){
	    for (vector<int>::iterator it = gene_table[i].begin(); it != gene_table[i].end(); ++it)
	      cout << " " << *it ;
	    cout << endl;
	  }
	  */
	  //derive permutation obtained from sorting survival_times
	  vector<int> permut(survival_times.size(), 0);
	  for (int i = 0 ; i != permut.size() ; i++) {
	    permut[i] = i;
	  }
	  sort(permut.begin(), permut.end(), compare_surv_times );
	  /*
	  for (int i = 0 ; i != permut.size() ; i++) {
	    cout << permut[i] << endl;
	  }
	  */
	  //now sort censoring vector using the permutation derived above
	  vector<int> censoring;
	  for (int i = 0; i < permut.size(); i++){
	    censoring.push_back( censoring_tmp[permut[i]] );
	  } 
	  //print content censoring
	  /*
	  cout << "Sorted censoring" << endl;
	  for (vector<int>::iterator it = censoring.begin(); it != censoring.end(); ++it)
	    cout << ' ' << *it;
	  cout << endl;
	  */
	  //find p_threshold: if a p-value estimated using the MC method is < p_threshold, then the exhaustive method or FPTAS are used
	  //to compute the p-value
	  double p_threshold;
	  string tmp = string( "-p_thres" );
	  int index_p_thres = -1;
	  for (int i=0; i<argc; i++){
	    if (tmp.compare(argv[i]) == 0){
	      index_p_thres = i;
	    }
	  }
	  if (index_p_thres > 0){
	     p_threshold = atof(argv[index_p_thres+1]);
	  }
	  else{
	    p_threshold = 1.0/double(number_iter);
	  }
	  cout << "[p_thresh: " << p_threshold << "]" << endl;
	  
	  //approximation factor
	  double e;
	  tmp = string( "-approx_factor" );
	  int index_e = -1;
	  for (int i=0; i<argc; i++){
	    if (tmp.compare(argv[i]) == 0){
	      index_e = i;
	    }
	  }
	  if (index_e > 0){
	     e = atof(argv[index_e+1]);
	  }
	  else{
	    e = 10.0;
	  }
	  cout << "[approximation factor: " << e << "]" <<  endl;
	  
	  //outfile name
	  ofstream outf;
	  tmp = string( "-out_file" );
	  int index_outf = -1;
	  for (int i=0; i<argc; i++){
	    if (tmp.compare(argv[i]) == 0){
	      index_outf = i;
	    }
	  }
	  if (index_outf > 0){
	     outf.open(argv[index_outf+1]);
	  }
	  else{
	    outf.open("output_ExaLT.txt");
	  }
	  
	  //write header on file
	  outf << "GENE_ID\tNUM_MUT_SAMPLE\tP_VALUE\tP_VALUE_LEFT\tP_VALUE_RIGHT\tMETHOD" << endl;
	  
	  //fix number of permutations for MC: depends on p_threshold, but make sure it is not too small
	  //for each geneID:
	  
	  double p_value;
	  double p_value_left;
	  double p_value_right;
	  
	  //tables to store info for later lookup
	  
	  //frequency[i] = mutation frequency for distributions in Vs_1, Prob_1, Vs_2, Prob_2
	  vector<int> frequency;
	  //vector of values for first side of FPTAS
	  vector<matrix_sim> Vs_1;
	  //vector of probabilities for first side of FPTAS
	  vector<matrix_sim> Prob_1;
	  //vector of values for second side of FPTAS
	  vector<matrix_sim> Vs_2;
	  //vector of probabilities for second side of FPTAS
	  vector<matrix_sim> Prob_2;
	  
	  for (int i=0; i<gene_table.size(); i++){
	    
	    //number of mutations for current gene
	    int freq_tmp = 0;
	    
	    //derive vector of mutations for the gene according to permutation above
	    vector<int> muts_tmp;
	    for (int j = 0; j < permut.size(); j++){
	      muts_tmp.push_back( gene_table[i][permut[j]] );
	      freq_tmp += gene_table[i][permut[j]];
	    } 
	    //print content censoring
	    /*
	    cout << "Sorted mutations for " << geneID[i] << "(freq=" << freq_tmp << ")" << endl;
	    for (vector<int>::iterator it = muts_tmp.begin(); it != muts_tmp.end(); ++it)
	      cout << ' ' << *it;
	    cout << endl;
	    */
	    //get p-value only if number of mutations is >= min_freq
	    if ( freq_tmp >= min_freq ){
	      //run MC to estimate p-value first
	      double p_value_MC;
	      double p_value_MC_left;
	      double p_value_MC_right;
	      //cout<<"Running quick estimate of the p-value..."<<endl;
	      string method = string("MC");
	      MC_estimate( muts_tmp, censoring, p_value_MC, p_value_MC_left, p_value_MC_right, number_iter );
	      if (p_value_MC == 0.0 ){
		      double min_pvalue_MC = 1.0/double(number_iter);
		      //cout <<"the estimated p-value is < "<<min_pvalue_MC<<endl;
		      p_value = min_pvalue_MC;
		      p_value_left = min_pvalue_MC;
 		      p_value_right = min_pvalue_MC;
	      }
	      else{
		      /*
		      cout <<"[The left p-value is estimated to be close to: "<<p_value_MC_left << "]" <<endl;
		      cout <<"[The right p-value is estimated to be close to: "<<p_value_MC_right << "]" <<endl;
		      cout<<"The p-value is estimated to be close to: "<<p_value_MC<<endl;
		      */
		      p_value = p_value_MC;
		      p_value_left = p_value_MC_right;
 		      p_value_right = p_value_MC_left;
	      }
	      //if p-value is < p_threshold:
	      if (p_value_MC <= p_threshold){
		//count number of combinations required by exact test
		/*
		cout <<"number of samples in group 0: "<<muts_tmp.size()-freq_tmp<<endl;
		cout <<"number of samples in group 1: "<<freq_tmp<<endl;
		*/
		unsigned long long ncomb = bcoeff(muts_tmp.size() , freq_tmp);
	  	//IF number combinations > threshold (fixed in code()): run FPTAS ELSE run exhaustive test
		if (ncomb <= 10000000){
		  //cout<<"Computing exact p-value..."<<endl;
		  method = string("exhaustive");
		  double p_value_exh;
		  double p_value_exh_left;
		  double p_value_exh_right;
		  exhaustive( muts_tmp, censoring, p_value_exh, p_value_exh_left, p_value_exh_right);
		  /*
		  cout <<"p-value (left) = " << p_value_exh_left << endl;
		  cout <<"p-value (right) = " << p_value_exh_right << endl;
		  cout<< geneID[i] << " p-value = "<<p_value_exh<<endl;
		  */
		  p_value = p_value_exh;
		  p_value_left = p_value_exh_left;
		  p_value_right = p_value_exh_right;
		}
		else{
		  method = string("FPTAS");
		  //check if the FPTAS has been run for this mutation frequency
		  int index_freq = -1;
		  for (int k=0; k<frequency.size(); k++){
		    if ( frequency[k] == freq_tmp ){
		      index_freq = k;
		    }
		  }
		  
		  if ( index_freq >= 0 ){
		    //cout << "****** USING LOOK-UP TABLES ******" << endl;
		    //compute the statistic
		    double V = 0;
		    double seen = 0;
		    for (int j=0;j<num_patients;j++) {
		      V = V + censoring[j] * (  (double)muts_tmp[j] -(freq_tmp-seen)/((double)num_patients-j));
		      seen += (double)muts_tmp[j];
		    }
		    //cout << "Statistic for first side: " << V << endl;
		    //now find the p-value
		    //first side
		    V=fabs(V);
		    double max_v_column_n1_index = -DBL_MAX;
		    double min_v_column_n1_index = DBL_MAX;
		    vector_max_min(Vs_1[index_freq].get_row(1), &max_v_column_n1_index, &min_v_column_n1_index);
		    double min_p_column_n1_index, max_p_column_n1_index;
 		    vector_max_min(Prob_1[index_freq].get_row(1), &max_p_column_n1_index, &min_p_column_n1_index);
		    /*
		    cout << "Prob_1[index_freq]: " << endl;
		    Prob_1[index_freq].display();
		    cout << "Vs_1[index_freq]: " << endl;
		    Vs_1[index_freq].display();
		    */
		    if ( V>max_v_column_n1_index ) {
		      p_value_right = 0.0;
		    }
		    else if ( V<=min_v_column_n1_index ) {
		      p_value_right = max_p_column_n1_index;
		    } else {
		      int l = intervalBound(Vs_1[index_freq], V, precision);
		      p_value_right=Prob_1[index_freq].get_element(1,l);
		    }
		    //cout << "p_value_right: " << p_value_right << endl;
		    //second side
		    V=-fabs(V);
		    max_v_column_n1_index = -DBL_MAX;
		    min_v_column_n1_index = DBL_MAX;
		    vector_max_min(Vs_2[index_freq].get_row(1), &max_v_column_n1_index, &min_v_column_n1_index);
		    vector_max_min(Prob_2[index_freq].get_row(1), &max_p_column_n1_index, &min_p_column_n1_index);
		    if ( V>max_v_column_n1_index ) {
		    //find the minimum p greater than 0
		      int size = Prob_2[index_freq].get_row(1).size();
		      //cout << "Size: " << size << endl;
		      double min_p = DBL_MAX;
		      for (int ii=0; ii<size; ++ii) {
			if ((Prob_2[index_freq].get_element(1,ii+1)<min_p) && (Prob_2[index_freq].get_element(1,ii+1)!=0)){
			  min_p = Prob_2[index_freq].get_element(1,ii+1);
			}
		      }
		      p_value_left = min_p;
		    }
		    else if ( V<=min_v_column_n1_index ) {
		      double e1 = 1.0 - pow( 1.0+e, -1.0/((double)num_patients) );
		      double tmp_dbl = (double)freq_tmp*log( (double)num_patients )/log(1+e1);
		      int k_max = tmp_dbl - (int)tmp_dbl > 0? (int)tmp_dbl+1 :  (int)tmp_dbl; //%maximum number of k
		      int k_max_index = k_max + 1;
		      p_value_left = Prob_2[index_freq].get_element(1,k_max_index);
		    } else {

		      int l = intervalBound_secondSide(Vs_2[index_freq], V, precision);
		      //l
		      p_value_left = Prob_2[index_freq].get_element(1,l);
		    }
		    //cout << "p_value_left: " << p_value_left << endl;
		    //total p_value
		    p_value = min(p_value_left+p_value_right, 1.0);
		    //cout << "p-value obtained from table: " << p_value << endl;
		  }
		  else {
		    frequency.push_back( freq_tmp );
    
		    matrix_sim v(1,1),p(1,1);
		    matrix_sim v2(1,1),p2(1,1);
		    matrix_sim v_tmp = v.get_col_as_row_matrix( freq_tmp+1 );
		    matrix_sim p_tmp = p.get_col_as_row_matrix( freq_tmp+1 );

		    double e2 = e;

		    pt_data data1;
		    data1.x = &muts_tmp;
		    data1.c = &censoring;
		    data1.e = &e;
		    data1.p_value = &p_value_right;
		    data1.v = &v;
		    data1.p = &p;

		    vector<int> censoring2 = censoring;
		    vector<int> muts_tmp2 = muts_tmp;

		    pt_data data2;
		    data2.x = &muts_tmp2;
		    data2.c = &censoring2;
		    data2.e = &e2;
		    data2.p_value = &p_value_left;
		    data2.v = &v2;
		    data2.p = &p2;
		    
		    //cout << "Obtaining controlled approximation..." << endl;
		    
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
		    v_tmp =  data1.v->get_col_as_row_matrix( freq_tmp+1 );
		    p_tmp =  data1.p->get_col_as_row_matrix( freq_tmp+1 );
		    Vs_1.push_back(v_tmp);
		    Prob_1.push_back(p_tmp);
		    v_tmp =  data2.v->get_col_as_row_matrix( freq_tmp+1 );
		    p_tmp =  data2.p->get_col_as_row_matrix( freq_tmp+1 );
		    Vs_2.push_back(v_tmp);
		    Prob_2.push_back(p_tmp);
		    p_value_left = *data2.p_value;
		    p_value_right = *data1.p_value;
		    p_value = min(p_value_left+p_value_right, 1.0);
		  }
		 }
	      }
	      //write output on file
	      outf << geneID[i] << "\t" << freq_tmp << "\t" << p_value << "\t" << p_value_left << "\t" << p_value_right << "\t" << method << endl;
	    }
	  }

	  outf.close();
	  
	  return 0;
	}



	//TODO: change output file names (to be passed as parameters)
	else{
	int line_c = 0;

        ifstream infile;
        infile.open(argv[1],ifstream::in);
        if(!infile.is_open()) {
                cerr<<"Error: Can not open data file " << argv[1] <<".\n" <<endl;
		cout<<endl;
                return -1;
        }
        string line, instring, cens_string;
        while ( getline(infile, line) ) {
                if (line_c == 0)
                        instring=line; //get 1st line
                else
                        cens_string=line; //get 2nd line, as sensor data
                line_c += 1;
        }
        infile.close();

        vector <int> x;
        int instring_len = instring.length();
        unsigned long long n_ones = 0;
        for (int i=0; i<instring_len; i++) {
                if (instring[i] == '0')
                        x.push_back(0);
                else {
                        n_ones += 1;
                        x.push_back(1);
                }
        }

        vector <int> c;
        int censtring_len = cens_string.length();
        for (int i=0; i<censtring_len; i++) {
                if (cens_string[i] == '0')
                        c.push_back(0);
                else {
                        c.push_back(1);
                }
        }

	vector<int> c2 = c;
	vector<int> x2 = x;
	double p_value;
	double p_value2;

	//run exhaustive if the number of combinations is not too large
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
		return 0;
	}

	//run MC to estimate p-value first
	double p_value_MC;
	double p_value_MC_left;
	double p_value_MC_right;
	cout<<"Running quick estimate of the p-value..."<<endl;
	MC_estimate( x, c, p_value_MC, p_value_MC_left, p_value_MC_right, number_iter );
	if (p_value_MC == 0.0 ){
		double min_pvalue_MC = 1.0/double(number_iter);
		cout <<"the estimated p-value is < "<<min_pvalue_MC<<endl;
	}
	else{
		cout <<"[The left p-value is estimated to be close to: "<<p_value_MC_left << "]" <<endl;
		cout <<"[The right p-value is estimated to be close to: "<<p_value_MC_right << "]" <<endl;
		cout<<"the p-value is estimated to be close to: "<<p_value_MC<<endl;
	}

	//now check if more controlled approximation is needed

	char input;
	do{
		cout << "Do you want a controlled approximation? [Y/N]" << endl;
		cin >> input;
	} while ( (input != 'Y') && (input != 'N') );
	if (input == 'N') {
                return 0;
	}
	else{

		cout << "Provide approximation factor (see README)"<<endl;
		double e;
		cin >> e;
		
		matrix_sim v(1,1),p(1,1);
		matrix_sim v2(1,1),p2(1,1);

		double e2 = e;

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
		data2.e = &e2;
		data2.p_value = &p_value2;
		data2.v = &v2;
		data2.p = &p2;

		cout << "Obtaining controlled approximation..." << endl;

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
		cout<<"p-value (left)  = "<<p_value<<endl;
		cout<<"p-value (right) = "<<p_value2<<endl;
		p_value = min(p_value+p_value2, 1);
		cout<<"p-value = "<<p_value<<endl;
		return 0;
	}
	}
}


