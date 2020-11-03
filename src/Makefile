# ROOTFLAG = `root-config --libs --cflags`
CFLAGS = -Wall -std=c++17 -ffast-math -O3
LIBFLAGS =  -I $(HOME)/nlopt/include -L$(HOME)/nlopt/lib -lnlopt -lm  -I /Applications/boost/boost_1_73_0 
all: 
	g++ $(CFLAGS) $(LIBFLAGS) -o ../bin/gfp_gaussian main.cpp

