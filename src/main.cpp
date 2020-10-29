#include "read_config.h"
#include "likelihood.h"
#include "utils.h"

#include <iostream> 
#include <iterator> 

int main(int argc, char** argv){

    Config config("config.txt");
    std::string infile = argv[1];
    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }
    // Read data
    std::vector<MOMAdata> cells =  getData(infile, 
                                        config.time_col,
                                        config.length_col,
                                        config.fp_col,
                                        config.delm);
    build_cell_genealogy(cells);

    std::cout << cells[0];

    Parameter_set random_params ;
    random_params.var_x = 1;

    // get the "trees" starting from all root cells
    std::vector<MOMAdata *> root_cells = get_roots(cells);

    for(long j=0; j<root_cells.size(); ++j){
        print_generation_tree(*root_cells[j], random_params);
        // double tl = total_likelihood(*root_cells[j], random_params);
    }
    
    std::cout << "Done." << std::endl;
    return 0;
}
