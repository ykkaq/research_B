
clean: 
	rm -f *.o


pdblas:
	g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 test_matrix.cpp -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp1
