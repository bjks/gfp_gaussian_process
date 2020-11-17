#include "CSVconfig.h"
#include "likelihood.h"
#include "minimizer_nlopt.h"
#include "tests.h"

#include <iostream> 
#include <iterator> 

int run_minimization(CSVconfig config, Parameter_set params, std::string infile){
    std::cout << params << "\n";

    /* Read data */
    std::cout << "-> Reading" << "\n";
    std::vector<MOMAdata> cells =  getData(infile, 
                                            config.time_col,
                                            config.length_col,
                                            config.fp_col,
                                            config.delm,
                                            config.cell_tags,
                                            config.parent_tags);

    /* genealogy */
    build_cell_genealogy(cells);
    // print_cells(cells);
    init_cells(cells, 5);

    std::cout << "-> Minimizaton" << "\n";
    /* minimization for tree starting from cells[0] */
    minimize_wrapper(&total_likelihood, cells[0], params);

    /* DEMO on how this works for all trees */
    /* get the "trees" starting from all root cells */
    std::vector<MOMAdata *> root_cells = get_roots(cells);
    std::vector<double> dummy_params;
    for(long j=0; j<root_cells.size(); ++j){
        print_generation_tree(dummy_params, *root_cells[j]);
    }
    return 0;
}


int main(int argc, char** argv){
    // test_likelihood();
    // return 0;

    CSVconfig config("csv_config.txt");
    Parameter_set params("parameter_bounds.txt");

    std::string infile = argv[1];    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }

    run_minimization(config, params, infile);
    std::cout << "Done." << std::endl;
    return 0;
}
