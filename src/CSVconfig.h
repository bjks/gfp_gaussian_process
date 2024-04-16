#include <fstream>
#include <iostream>
#include "utils.h"


class CSVconfig {
    /*  
    * Class for handling the configuration file passed on to the constructor 
    * that specifies the columns of the csv file read into MOMAdata instances as well as the delimiter
    */
public:
    // Defaults
    std::string time_col = "time";
    double rescale_time = 1.;
    std::string length_col = "length";
    bool length_islog = false;
    std::string fp_col = "gfp";
    double fp_auto = 0;
    std::string delm = ",";
    std::string segment_col = "";
    std::string filter_col = "";

    std::vector<std::string> cell_tags {"cell_id"};
    std::vector<std::string> parent_tags {"parent_id"};

    CSVconfig(std::string filename) {
        std::ifstream fin(filename);
        // fin.open(filename);
        std::string line;
        std::vector<std::string> parts;

        // Overwrite defaults if in config file
        while (getline(fin, line)) {
            if (line[0] != '#' && line.size()){
                parts = split_string_at(line, "=");

                // remove whitespaces from the ends
                parts[0] = trim(parts[0]);
                parts[1] = trim(parts[1]);

                if (parts[0] == "time_col"){
                    time_col = parts[1];
                }
                else if (parts[0] == "rescale_time"){
                    try{
                        rescale_time = std::stod(parts[1]);
                    }
                    catch(std::exception &e){
                        _file_log << "(CSVconfig) ERROR: rescale_time in 'csv_file' cannnot be processed (" << e.what() <<")" << std::endl;
                        throw;
                    }
                }
                else if (parts[0] == "length_col"){
                    length_col = parts[1];
                }
                else if (parts[0] == "length_islog"){
                    length_islog = string2bool(parts[1]);
                }
                else if (parts[0] == "fp_col"){
                    fp_col = parts[1];
                }
                else if (parts[0] == "fp_auto"){
                    try{
                        fp_auto = std::stod(parts[1]);
                    }
                    catch(std::exception &e){
                        _file_log << "(CSVconfig) ERROR: fp_auto in 'csv_file' cannnot be processed (" << e.what() <<")" << std::endl;
                        throw;
                    }
                }
                else if (parts[0] == "delm"){
                    delm = parts[1];
                }
                else if (parts[0] == "cell_tags"){
                    cell_tags.clear();
                    std::vector<std::string> val_split;
                    val_split = split_string_at(parts[1], ",");
                    for (size_t i=0; i<val_split.size(); ++i){
                        cell_tags.push_back(trim(val_split[i]));
                    }
                }
                else if (parts[0] == "parent_tags"){
                    parent_tags.clear();
                    std::vector<std::string> val_split;
                    val_split = split_string_at(parts[1], ",");
                    for (size_t i=0; i<val_split.size(); ++i){
                        parent_tags.push_back(trim(val_split[i]));
                    }
                }
                else if (parts[0] == "segment_col"){
                    segment_col = parts[1];
                }
                else if (parts[0] == "filter_col"){
                    filter_col = parts[1];
                }
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const CSVconfig& config);

};

std::ostream& operator<<(std::ostream& os, const CSVconfig& config){
    int col = 15;
    os << "Configuration used for reading the input file\n"; 
    os << "_____________________________________________\n"; 
    os          << pad_str("time_col:", col)  << config.time_col << "\n" 
                << pad_str("rescale_time:", col)  << config.rescale_time << "\n" 
                << pad_str("length_col:", col) <<  config.length_col << "\n" 
                << pad_str("length_islog:", col) <<  config.length_islog << "\n" 
                << pad_str("fp_col:", col) << config.fp_col << "\n" 
                << pad_str("fp_auto:", col)  << config.fp_auto << "\n" 
                << pad_str("delm:", col) << config.delm << "\n"; 

    if (!config.segment_col.empty())
        os << pad_str("segment_col:", col) << config.segment_col << "\n";
    if (!config.filter_col.empty())
        os << pad_str("filter_col:", col) << config.filter_col << "\n";
    
    os << pad_str("cell_tags:", col) ;
    for(size_t i=0; i < config.cell_tags.size(); i++){
        os << config.cell_tags[i] << ' ';
    }
    os << "\n"; 
    os << pad_str("parent_tags:", col);
    for(size_t i=0; i < config.parent_tags.size(); i++){
        os << config.parent_tags[i] << ' ';
    }
    os << "\n";
    return os;
}