#include "CSVconfig.h"
#include "likelihood.h"
#include "utils.h"
#include "maximize_likelihood.h"

#include <iostream> 
#include <iterator> 

int main(int argc, char** argv){

    CSVconfig config("csv_config.txt");
    Parameter_set params("parameter_bounds.txt");

    std::cout << params << "\n";

    numerical_minimization(params);

    // Read data
    std::string infile = argv[1];    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }
    std::vector<MOMAdata> cells =  getData(infile, 
                                            config.time_col,
                                            config.length_col,
                                            config.fp_col,
                                            config.delm);
    build_cell_genealogy(cells);


    // get the "trees" starting from all root cells
    std::vector<MOMAdata *> root_cells = get_roots(cells);

    for(long j=0; j<root_cells.size(); ++j){
        print_generation_tree(*root_cells[j], params);
    }
    
    std::cout << "Done." << std::endl;
    return 0;
}
