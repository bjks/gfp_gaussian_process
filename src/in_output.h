#include <fstream>
#include <iostream>
#include <sstream>  

#include "Parameters.h"

std::string _infile;
std::string _outfile;
int _iteration = 0;
int _print_level = 0;


std::string outfile_name_minimization(std::string infile, Parameter_set params){
    std::vector<std::string> infile_split, base_split;
    infile_split = split_string_at(infile, "/");
    base_split = split_string_at(infile_split[infile_split.size()-1], ".");

    std::string outfile = "";
    for(int i=0; i<infile_split.size()-1; ++i){
        outfile += infile_split[i] + "/";
    }
    outfile += base_split[0] + "_out/";
    std::__fs::filesystem::create_directory(outfile);
    outfile += base_split[0]+ "_f";
    for(int i=0;i<params.all.size();++i){
        if (!params.all[i].bound && !params.all[i].fixed){
            outfile += std::to_string(i);
        }
    }

    outfile += "_b";
    for(int i=0;i<params.all.size();++i){
        if (params.all[i].bound){
            outfile += std::to_string(i);
        }
    }

    return outfile + ".csv";
}

std::string outfile_name_scan(std::string infile, std::string var){
    std::vector<std::string> infile_split, base_split;
    infile_split = split_string_at(infile, "/");
    base_split = split_string_at(infile_split[infile_split.size()-1], ".");

    std::string outfile = "";
    for(int i=0; i<infile_split.size()-1; ++i){
        outfile += infile_split[i] + "/";
    }
    outfile += base_split[0] + "_out/";
    std::__fs::filesystem::create_directory(outfile);
    outfile += base_split[0]+ "_scan_" + var;
    return outfile + ".csv";
}

std::string outfile_name_prediction(std::string infile){
    std::vector<std::string> infile_split, base_split;
    infile_split = split_string_at(infile, "/");
    base_split = split_string_at(infile_split[infile_split.size()-1], ".");

    std::string outfile = "";
    for(int i=0; i<infile_split.size()-1; ++i){
        outfile += infile_split[i] + "/";
    }
    outfile += base_split[0] + "_out/";
    std::__fs::filesystem::create_directory(outfile);
    outfile += base_split[0]+ "_prediction";
    return outfile + ".csv";
}

void setup_outfile_params(std::string outfile, Parameter_set params){
    params.to_csv(outfile);
    std::ofstream file(outfile,std::ios_base::app);

    file << "\nlikelihoods:\niteration,";
    for (int i=0; i<params.all.size(); ++i){
        file << params.all[i].name << ",";
    }
    file << "likelihood" <<"\n";
    file.close();
}

void setup_outfile_mean(std::string outfile, Parameter_set params){
    params.to_csv(outfile);
    std::ofstream file(outfile,std::ios_base::app);

    file << "\nmean_x,mean_g,mean_lambda,mean_q\n";
    file.close();
}

