#include "CSVconfig.h"
#include "likelihood.h"
#include "utils.h"
#include "minimizer_nlopt.h"

#include <iostream> 
#include <iterator> 

int main(int argc, char** argv){

    CSVconfig config("csv_config.txt");
    Parameter_set params("parameter_bounds.txt");

    std::cout << params << "\n";

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
    minimize_wrapper(&total_likelihood, cells[0], params);

    print_cells(cells);

    // get the "trees" starting from all root cells
    std::vector<MOMAdata *> root_cells = get_roots(cells);
    std::vector<double> dummy_params;
    for(long j=0; j<root_cells.size(); ++j){
        print_generation_tree(dummy_params, *root_cells[j]);
    }
    
    std::cout << "Done." << std::endl;
    return 0;
}
