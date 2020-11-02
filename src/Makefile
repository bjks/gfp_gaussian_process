CFLAGS= -Wall -I /Applications/boost/boost_1_73_0 -std=c++17 -ffast-math -O3
ROOTFLAG = `root-config --libs --cflags`

all: 
	g++ $(CFLAGS) $(ROOTFLAG) -o ../bin/gfp_gaussian main.cpp

