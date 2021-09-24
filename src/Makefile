# CFLAGS = -std=c++17 -ffast-math -O3
CFLAGS = -std=c++17 -O3

FADFLAG =  Faddeeva.cc

NLOPTLAGS = -I ~/nlopt/include -L ~/nlopt/lib -lnlopt -lm  
NLOPTLAGS64 = -I ~/nlopt/include -L ~/nlopt/lib64 -lnlopt -lm  

BOOSTFLAG = -I /Applications/boost/boost_1_73_0 -L /opt/local/lib  -lboost_iostreams-mt -lz

EIGENFLAG = -I ~/eigen  

local: 
	g++ main.cpp -Wall $(CFLAGS) $(NLOPTLAGS) $(BOOSTFLAG) $(EIGENFLAG) $(FADFLAG) -o ../bin/gfp_gaussian

variation:
	g++ main.cpp -Wall $(CFLAGS) $(NLOPTLAGS) $(EIGENFLAG) $(FADFLAG) -o ../bin/gfp_gaussian_variation

cluster:
	ml Eigen/3.3.7; ml Boost/1.74.0-foss-2018b; g++ main.cpp $(CFLAGS) $(NLOPTLAGS64) $(FADFLAG) -static-libstdc++ -lstdc++fs -o ../bin/gfp_gaussian







