# ROOTFLAG = `root-config --libs --cflags`
CFLAGS = -Wall -std=c++17 -ffast-math -O3
NLOPTLAGS = -I $(HOME)/nlopt/include -L$(HOME)/nlopt/lib -lnlopt -lm  
EIGENFLAG = -I /Applications/Eigen/eigen  

all: 
	g++ $(CFLAGS) $(NLOPTLAGS) $(BOOSTFLAG) $(EIGENFLAG) -o ../bin/gfp_gaussian main.cpp

