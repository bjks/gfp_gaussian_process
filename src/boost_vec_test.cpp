#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


void append_vec(boost::numeric::ublas::vector<double> &v, double elem){
    v.resize(v.size()+1);
    v[v.size()-1] = elem;
}

int main(int argc, char** argv){
    
    boost::numeric::ublas::vector<double> vec;

    append_vec(vec, 2);
    append_vec(vec, 3);
    append_vec(vec, 4);

    for (double v: vec){
        std::cout << v << std::endl;
    }

    std::cout << "Done." << std::endl;
    return 0;
}
