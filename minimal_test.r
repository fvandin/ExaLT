#requires R package Rcpp, calling compiled cpp code from R
library(Rcpp)

#load the compiled cpp dynamic link library
dyn.load('R_ExaLT.so');

#call the function, named 'R_FPTAS', with 2 arguments, 1. table data file name, 2, column
.Call('R_ExaLT', 'minimal_example_R.txt', 'HIPK2')
