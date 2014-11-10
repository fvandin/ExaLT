ExaLT
=====

Exact Log-rank Test - Accurate Computation of Survival Statistics in Genome-wide Studies

Copyright 2010,2011,2012,2013,2014 Fabio Vandin, Brown University, Providence, RI.

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


ExaLT - Exact Log-rank Test

Version: 1.1

===========================================================

COMPILING AND RUNNING ExaLT (C++ VERSION)

First try to compile the multithreaded version of ExaLT 
use the following commands:
$ cd /directory/with/ExaLT/code/
$ make main_FPTAS

The multithreaded version of ExaLT requires pthread.

If that does not produce the file 'ExaLT', or if you prefer to use the single
threaded version, use the following commands:
$ cd /directory/with/ExaLT/code/
$ make main_FPTAS_noMultiThread

- Running ExaLT

ExaLT can be run in two modes: in the first one, it will test the difference in
survival between two populations; in the second one, it will test the association
between survival time and a number of (binary) features provided as columns of a
table (each feature is tested separately). See the two sections below for details.

- Running ExaLT: single test

To run ExaLT for a single test, use:
./ExaLT input_file

where input_file is the file describing the order of events in the two populations
and the censoring information (see "C++ INPUT FORMAT - Single test" below for
details on the input format).
The output is printed on screen.

ExaLT uses the exhaustive computation if the sample space is small enough that
the running time is reasonable. For example, running:

./ExaLT minimal_example.txt 

will print on screen:

number of samples in group 0: 3
number of samples in group 1: 1
Computing exact p-value...
p-value (left)   = 0.75
p-value (right)  = 0.25
p-value = 1

"p-value (right)" is the pvalue for the right part of the tail; "p-value (left)"
is the pvalue for the left part of the tail; "p-value" is the two-tail p-value.

When the sample space is too large, ExaLT first computes an estimate of the
p-value using a Monte Carlo approach. While pretty accurate, this estimate may be
non conservative; then ExaLT asks if a conservative, controlled estimated is
required. For example, it will print some text like the following:
	
Running quick estimate of the p-value...
[The left p-value is estimated to be close to: 0.053]
[The right p-value is estimated to be close to: 0.049]
The p-value is estimated to be close to: 0.12
Do you want a controlled approximation? [Y/N]

If Y is chosen, then ExaLT ask to provide an approximation factor (must be > 0).
If the approximation factor provided is E and the exact p-value is P, then ExaLT
returns a value P* such that P <= P* < (1+E)P. That is, P* is at most a factor
(1+E) of the exact p-value. The same guarantees are provided for the left and
the right p-value. We suggest to set E=10 if in doubt, since the computation may
require longer time for smaller values of E.

Do you want a controlled approximation? [Y/N]
Y
Provide approximation factor (see README)
10      
Obtaining controlled approximation...
p-value (left)  = 0.0532
p-value (right) = 0.0482
p-value = 0.114

(Note that the values above are not obtained with the small
minimal_example.txt!)

- Running ExaLT: table

In the second mode, ExaLT will run on a table (see "C++ INPUT FORMAT - Table" section
below). ExaLT uses the genes (or any binary feature) in the table to split the
patients/samples into two populations, and for each split computes the p-value of the
difference among the survival distribute of the two populations.

To run ExaLT in this mode, use:

./ExaLT table_file_name -table min_frequency [-p_thres VALUE -out_file NAME -approx_factor VALUE]

where required input parameters are:
- table_file_name: name of file with input table (see "C++ INPUT FORMAT - Table"
section);
- -table: required argument (specifies that table_file_name contains a table);
- min_frequency: minimum frequency for of mutation of a gene to be included among the
tested genes.

Optional parameters are:
- -p_thres VALUE: threshold for p-values; tests with MC estimated p-values below this
thresholds are then run using the exhaustive test or the controlled approximation
method (default: 0.00001)
-out_file NAME: name of the output file (default: output_ExaLT.txt)
-approx_factor VALUE: approximation factor for the controlled approximation
(default: 10)

The output format is given by the results of the tests with one gene per line,
as follows:
GENE_ID NUM_MUT_SAMPLE  P_VALUE P_VALUE_LEFT    P_VALUE_RIGHT   METHOD

where:
- GENE_ID: name of the gene (obtained from table_file_name);
- NUM_MUT_SAMPLE: mutation frequency of the gene (i.e., number of samples with mutations
in the gene);
- P_VALUE: p-value for the association between mutation status of the gene and survival
time;
- P_VALUE_LEFT: left p-value for the association between mutation status of the gene and
survival time;
- P_VALUE_RIGHT: right p-value for the association between mutation status of the gene
and survival time;
- METHOD: method that has been used to compute the p-value: exhaustive = exhaustive
algorithm; MC = Monte Carlo estimate; FPTAS = controlled approximation algorithm.

See section "EXAMPLE OF AN INPUT FILE - Table" for an example.

NOTE: ties in survival times are handled by using an arbitrary order.

===========================================================

COMPILING AND RUNNING THE R VERSION

These are the instructions to compile ExaLT code and install the necessary
software to call it from R.

NOTE: 1 and 2 are to get Rcpp package installed. If the corresponding versions R
and Rcpp are installed in a different way are aworking together, that is fine,
and you can skip to step 3.

1. Need R version >2.15, to use Rcpp package.
Download most recent version of R:
http://cran.r-project.org/bin/linux/ubuntu/README

2. Download Rcpp (to interface between cpp and R).
Instructions to install an R package:
http://math.usask.ca/~longhai/software/installrpkg.html
Let /my/own/Rccp_dir be the folder with the Rccp downloaded files.
$ R CMD INSTALL Rcpp -l /my/own/Rcpp_dir/
(On some systems, it needs 'sudo', because it will install the package in
/usr/lib/...)

3. For  "R CMD SHLIB conv.cpp" to work correctly
In ~/.bashrc (or corresponding file, e.g. ~/.bash_profile on OS X) add:
export PKG_CPPFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`
export PKG_CPPFLAGS=" -I. -lpthread "$PKG_CPPFLAGS
export PKG_LIBS=`Rscript -e "Rcpp:::LdFlags()"`

4. Compile the C++ code to be called by R:
$ cd /directory/with/ExaLT/code/
$ make

5. Test the code:
$ R
> source('test.r')

To call ExaLT from R, use the following in R (see also test.r):

#load R package Rcpp, calling compiled cpp code from R
>library(Rcpp)

#load the compiled cpp dynamic link library; the following assumes you are in
the dir with the C++ compiled code
>dyn.load('R_ExaLT.so');

#call the function, named 'R_ExaLT', with 2 arguments: 1. table data file name,
2. column with feature defining the two groups (with values 0 and 1). In the
example, 'HIPK2' is the name of the feature.
.Call('R_ExaLT', 'tableR2.txt', 'HIPK2')

Note that running ExaLT from R will correspond to run the C++ version of 'Single text'
on each gene separately (see "COMPILING AND RUNNING ExaLT (C++ VERSION)" section for 
details on execution and output).

===========================================================

C++ INPUT FORMAT

- Single test

The input file for C++ should contain two lines that specify the survival data
(time, censoring) for the sorted order of patients:
1) the first line specifies a sequence reporting the population (0 or 1) in
which the event (censored or not is found). For example, if the first event is
found in population 0, the second in population 0, the third in population 1,
and the fourth in population 0, the first line is:
0010
2) the second line specifies the censoring parameter: again a sequence of 0's
and 1's, with 0 in the i-th position if the i-th event is censored, and 1
otherwise. For example, if the second and third patients (after ordering by
increasing survival time) are censored the second line is:
1001

That is, in this example the input is given by:
0010
1001

- Table

The input table uses rows to define patients, and columns to define genes (or any binary
feature). The first row is used to define the gene names and has the following format
(tokens are TAB separated):
patientID	survival_time   censoring_status  GENE1 GENE2 GENE3 ?
in this case, ExaLT will consider GENE1, GENE2, and GENE3 as gene names. (Note the
patientID, survival_time, and censoring_time can be substituted by any other string.)
From the second row to the last one, the format is:
patient_for_row	survival_time   censoring_status  MUT_IN_GENE1? MUT_IN_GENE2? MUT_IN_GENE3? ?
where MUT_IN_GENE1? is 0 if GENE1 is not mutated in the patient described by the row,
and 1 if GENE1 is mutated in the patient described by the row, and analogously for the
other genes.

For an example, see "EXAMPLE OF INPUT FILE - Table".

===========================================================

R INPUT FORMAT

* tableR.txt, tableR2.txt, and tableR3.txt are other small, example input tables
to try. The input format for a R table to be used with ExaLT is:
- header: patientID\tSTATUStTIME\tGENE_NAME1\tGENE_NAME2\t...
Note that instead of GENE_NAME1,GENE_NAME2,... general features names can be
used; the feature is used only to identify 2 groups: one with value 0 for the
feature (e.g., GENE_NAME1 is 0 if gene is not mutated in the patient); one with
value 1 for the feature (e.g., GENE_NAME1 is 1 if gene is mutated in the
patient). To have the code to work properly, the group with value 0 for the
feature should be the largest group (i.e., the number of lines with value 0 for
the feature should be larger than the number of lines with value 1 for the
feature).
- for each row: patientID\tCensoringStatus\tSurvivalTime\tValue1\tValue2...
The first value is the patientID; the second is the status (for censoring),
either DECEASED or LIVING (1 and 0 can also be used in alternative, respectively
for DECEASED and LIVING); the third is the survival time; from the forth, values
of features (in order) are reported. Values are TAB (\t) separated. See also
example below ("INPUT FORMAT")

===========================================================

EXAMPLE OF INPUT FILE

- Single test (R and C++)

The following is a minimal example of the input for R. In the example there are
4 patients, and the only feature is the mutation status of gene HIPK2 (1 =
mutated, 0 = not mutated).

PATIENTID	STATUS	TIME	HIPK2
sample1	LIVING	90.5	0
sample2	LIVING	108.0	1
sample3	DECEASED	150.1	0
sample4	DECEASED	87.0

See also minimal_example_R.txt

The corresponding input file for the C++ version is obtained by sorting the
patients by increasing survival time, and then spelling out the value of the
HIPK2 feature on the first line, and the STATUS on the second line after
substituting LIVING with 0 and DECEASED with 1:

0010
1001

See also minimal_example.txt

- Table

see table_example.txt in the folder containing the source code of ExaLT. The content of
the file is:

patientID       survival_time   censoring_status  TP53    PIK3CA  PTEN
pa1     5.2     1       1       0       1
ce1     4.2     0       0       0       0
pa2     6.2     1       1       0       1
pa3     7.2     1       0       0       0
ce2     6.4     0       0       1       0

In this example there are 5 patients (IDs:pa1,ce1,pa2,pa3, and ce2) and 3 genes (TP53,
PIK3CA, and PTEN). Patient pa1 has survival time 5.2, and IS NOT censored (censoring_status
is 1). Patient ce1 has survival time 4.2, and IS censored (censoring_status is 0). Gene
TP53 is mutated in patients pa1 and pa2, while PIK3CA is mutated in patient ce2. 

Running the following command (contained in test_table.sh) in the folder containing the
source code of ExaLT:

./ExaLT table_example.txt -table 1 -p_thres 0.99 -out_file output_test.txt -approx_factor 1

executes the test on all genes mutated in >= 1 patients (-table 1), uses the exhaustive or
controlled approximation for every gene with an MC p-value estimated p-value < 0.99, uses
and approximation factor of 1 (-approx_factor 1), and writes results to file output_test
(-out_file output_test.txt).

After running the command above, the content of output_test.txt is:

GENE_ID NUM_MUT_SAMPLE  P_VALUE P_VALUE_LEFT    P_VALUE_RIGHT   METHOD
TP53    2       0.1     0       0.1     exhaustive
PIK3CA  1       0.6     0.4     0.2     exhaustive
PTEN    2       0.1     0       0.1     exhaustive

For example, the second row tells that TP53, mutated in 2 samples, has p-value 0.1 (left
p-value = 0, right p-value = 0.1). In this very small example, all p-values have been
computed by the exhaustive method.

===========================================================

REMOVE COMPILED FILES

To remove the compiled files, just run 
$ cd /directory/with/ExaLT/code/
$ make clean

Note that the above does not remove the source (just delete the folder
/directory/with/ExaLT/code/ to remove the sources).

REFERENCES:

If you use ExaLT in your research, please cite:
F. Vandin, A. Papoutsaki, B.J. Raphael, and E. Upfal.
Genome-Wide Survival Analysis of Somatic Mutations in Cancer.
In Proceedings 17th International Conference on Research in Computational Molecular Biology
(RECOMB), LNCS 7821, 2013, pp 285-286. 

CONTACT:

Fabio Vandin: vandinfa@cs.brown.edu