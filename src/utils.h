#include <iostream>
#include <vector>
#include <string>

#ifndef UTILS_H
#define UTILS_H

void pvector(std::vector <std::string> const &a) {
    for(int i=0; i < a.size(); i++){
        std::cout << a[i] << ' ';
    }
    std::cout << std::endl;
}

void pvector(std::vector <double> const &a) {
    for(int i=0; i < a.size(); i++){
        std::cout << a[i] << ' ';
    }
    std::cout << std::endl;
}

std::vector<std::string>  split_string_at(std::string s, std::string delimiter=","){
    int pos = 0;
    std::vector<std::string> splitted;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        splitted.push_back(s.substr(0, pos));
        s.erase(0, pos + delimiter.length());
    }
    splitted.push_back(s);
    return splitted;
}

std::string trim(const std::string& str, char trim_char =' '){
    int first = str.find_first_not_of(trim_char);
    if (std::string::npos == first){
        return str;
    }
    int last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

#endif 
