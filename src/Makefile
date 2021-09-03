# CFLAGS = -std=c++17 -ffast-math -O3
CFLAGS = -std=c++17 -O3

FADFLAG =  Faddeeva.cc

NLOPTLAGS = -I ~/nlopt/include -L ~/nlopt/lib -lnlopt -lm  
NLOPTLAGS64 = -I ~/nlopt/include -L ~/nlopt/lib64 -lnlopt -lm  

EIGENFLAG = -I ~/eigen  

local: 
	g++ main.cpp -Wall $(CFLAGS) $(NLOPTLAGS) $(EIGENFLAG) $(FADFLAG) -o ../bin/gfp_gaussian

# variation:
# 	g++ main.cpp -Wall $(CFLAGS) $(NLOPTLAGS) $(BOOSTFLAG) $(EIGENFLAG) $(FADFLAG) -o ../bin/gfp_gaussian_variation

cluster:
	ml GCC/8.3.0; ml Eigen/3.3.7; g++ main.cpp $(CFLAGS) $(NLOPTLAGS64)  $(EIGENFLAG) $(FADFLAG) -static-libstdc++ -lstdc++fs -o ../bin/gfp_gaussian







