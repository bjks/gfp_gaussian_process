#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>


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

    CSVconfig(std::string filename) {
        std::ifstream fin(filename);
        // fin.open(filename);
        std::string line;
        std::vector<std::string> vec;

        if(fin) {
            // Overwrite defaults if in config file
            while (getline(fin, line)) {
                if (line[0] != '#' && line.size()){
                    boost::algorithm::split(vec, line, boost::is_any_of("="));

                    // remove whitespaces from the ends
                    boost::algorithm::trim(vec[0]);
                    boost::algorithm::trim(vec[1]);

                    if (vec[0] == "time_col"){
                        time_col = vec[1];
                    }
                    else if (vec[0] == "length_col"){
                        length_col = vec[1];
                    }
                    else if (vec[0] == "fp_col"){
                        fp_col = vec[1];
                    }
                    else if (vec[0] == "delm"){
                        delm = vec[1];
                    }
                }
            }
        }
        else {
            std::cout<<"Config file NOT found, use defaults" << std::endl; 
        }

        std::cout << "--- Configuration used for reading the input file: ---" << std::endl; 
        std::cout   << "# time_col: "  << time_col << "\n" \
                    << "# length_col: " <<  length_col << "\n" \
                    << "# fp_col: " << fp_col << "\n" \
                    << "# delm: " << delm << "\n" ;
        std::cout << "------------------------------------------------------" << std::endl; 
        
    }
};


