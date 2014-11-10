test_FPTAS: R_ExaLT.cpp FPTAS.cpp logrank_exact.cpp MC_logrank.cpp  matrix_sim.cpp  table.cpp FPTAS.h  logrank_exact.h MC_logrank.h  matrix_sim.h  table.h
	R CMD SHLIB R_ExaLT.cpp FPTAS.cpp logrank_exact.cpp MC_logrank.cpp  matrix_sim.cpp  table.cpp
main_FPTAS: 
	g++ -O3 -c main_FPTAS.cpp logrank_exact.cpp MC_logrank.cpp  matrix_sim.cpp  table.cpp
	g++ -o ExaLT -lpthread main_FPTAS.o logrank_exact.o MC_logrank.o  matrix_sim.o  table.o 
main_FPTAS_noMultiThread:
	g++ -O3 -c main_FPTAS_noMultiThread.cpp logrank_exact.cpp MC_logrank.cpp  matrix_sim.cpp  table.cpp
	g++ -o ExaLT main_FPTAS_noMultiThread.o logrank_exact.o MC_logrank.o  matrix_sim.o  table.o
clean:
	rm R_ExaLT.so table.o matrix_sim.o MC_logrank.o logrank_exact.o main_FPTAS_noMultiThread.o main_FPTAS.o FPTAS.o R_ExaLT.o ExaLT
