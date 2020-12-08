CFLAGS = -Wall -std=c++14 -ffast-math -O3
FADFLAG =  Faddeeva.cc

NLOPTLAGS = -I $(HOME)/nlopt/include -L$(HOME)/nlopt/lib -lnlopt -lm  
# NLOPTLAGS = -I nlopt/include -L nlopt/lib -lnlopt -lm  
EIGENFLAG = -I $(HOME)/eigen  

compile_all: 
	g++ $(CFLAGS) $(NLOPTLAGS) $(BOOSTFLAG) $(EIGENFLAG) $(FADFLAG) -o ../bin/gfp_gaussian main.cpp

