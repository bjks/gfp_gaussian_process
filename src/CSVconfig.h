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
    std::string time_col = "time_sec";
    std::string length_col = "length_um";
    std::string fp_col = "gfp_nb";
    std::string delm = ",";

    std::vector<std::string> cell_tags {"date", "pos", "gl", "id"};
    std::vector<std::string> parent_tags {"date", "pos", "gl", "parent_id"};

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
                else if (parts[0] == "length_col"){
                    length_col = parts[1];
                }
                else if (parts[0] == "fp_col"){
                    fp_col = parts[1];
                }
                else if (parts[0] == "delm"){
                    delm = parts[1];
                }
                else if (parts[0] == "cell_tags"){
                    cell_tags.clear();
                    std::vector<std::string> val_split;
                    val_split = split_string_at(parts[1], ",");
                    for (int i=0; i<val_split.size(); ++i){
                        cell_tags.push_back(trim(val_split[i]));
                    }
                }
                else if (parts[0] == "parent_tags"){
                    parent_tags.clear();
                    std::vector<std::string> val_split;
                    val_split = split_string_at(parts[1], ",");
                    for (int i=0; i<val_split.size(); ++i){
                        parent_tags.push_back(trim(val_split[i]));
                    }
                }
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const CSVconfig& config);

};

std::ostream& operator<<(std::ostream& os, const CSVconfig& config){
    std::cout << "Configuration used for reading the input file\n"; 
    std::cout << "_____________________________________________\n"; 
    std::cout   << "\ttime_col: "  << config.time_col << "\n" \
                << "\tlength_col: " <<  config.length_col << "\n" \
                << "\tfp_col: " << config.fp_col << "\n" \
                << "\tdelm: " << config.delm << "\n" \
                << "\tcell_tags: " ;
    pvector(config.cell_tags);
    std::cout << "\tparent_tags: " ;
    pvector(config.parent_tags);
    std::cout << "\n\n"; 
    return os;
}