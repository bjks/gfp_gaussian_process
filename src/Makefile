CFLAGS = -Wall -std=c++17 -ffast-math -O3
FADFLAG =  Faddeeva.cc

NLOPTLAGS = -I ~/nlopt/include -L ~/nlopt/lib -lnlopt -lm  
NLOPTLAGS64 = -I ~/nlopt/include -L ~/nlopt/lib64 -lnlopt -lm  

EIGENFLAG = -I ~/eigen  

local: 
	g++ main.cpp $(CFLAGS) $(NLOPTLAGS) $(BOOSTFLAG) $(EIGENFLAG) $(FADFLAG) -o ../bin/gfp_gaussian

cluster:
	g++ main.cpp $(CFLAGS) $(NLOPTLAGS64) $(BOOSTFLAG) $(EIGENFLAG) $(FADFLAG) -lstdc++fs -o ../bin/gfp_gaussian







