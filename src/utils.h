#include <iostream>
#include <vector>
#include <string>
#include <sstream>  
#include <map> 
#include <filesystem>


#ifndef UTILS_H
#define UTILS_H

std::string pad_str(std::string s, const size_t num, const char paddingChar = ' '){
    /* pads string with paddingChar to match desired length */
    if(num > s.size())
        s.insert(s.end(), num - s.size(), paddingChar);
    return s;
}

std::string pad_str(double d, const size_t num, const char paddingChar = ' '){
    /* pads string conversion of double with paddingChar to match desired length */
    std::stringstream buffer;
    buffer << d;
    std::string s = buffer.str();
    if(num > s.size())
        s.insert(s.end(), num - s.size(), paddingChar);
    return s;
}

void pvector(std::vector <std::string> const &a) {
    for(size_t i=0; i < a.size(); i++){
        std::cout << a[i] << ' ';
    }
    std::cout << std::endl;
}

void pvector(std::vector <double> const &a) {
    for(size_t i=0; i < a.size(); i++){
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
    /* trims string from whitespaces */
    int first = str.find_first_not_of(trim_char);
    if (std::string::npos == first){
        return str;
    }
    int last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    /* numpy like arange */
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}


std::string default_out_dir(std::string infile){
    /* composes the default out_dir */
    std::vector<std::string> infile_split, base_split;
    infile_split = split_string_at(infile, "/");
    base_split = split_string_at(infile_split[infile_split.size()-1], ".");

    std::string out_dir = "";
    for(size_t i=0; i<infile_split.size()-1; ++i){
        out_dir += infile_split[i] + "/";
    }
    out_dir += base_split[0] + "_out/";
    return out_dir;
}

std::string out_dir(std::map<std::string, std::string> arguments){
    /* returns the given out_dir name (ending with '/') or returns the default out_dir */
    std::string out_dir;
    if (!arguments.count("outdir"))
        out_dir = default_out_dir(arguments["infile"]);
    else{
        out_dir = arguments["outdir"];
        if (out_dir.back() != '/')
            out_dir += "/";
    }
    std::filesystem::create_directory(out_dir);
    return out_dir;
}


std::string file_base(std::string infile){
    /* returns the base of the file name (without dir name) */
    std::vector<std::string> infile_split, base_split;    
    infile_split = split_string_at(infile, "/");
    base_split = split_string_at(infile_split[infile_split.size()-1], ".");
    return base_split[0];
}

#endif 
