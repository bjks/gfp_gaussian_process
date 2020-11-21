#include "CSVconfig.h"
#include "likelihood.h"
#include "minimizer_nlopt.h"
#include "tests.h"

#include <iostream> 
#include <iterator> 


void run_minimization(CSVconfig config, Parameter_set params, std::string infile){
    /* Read data */
    std::cout << "-> Reading" << "\n";
    std::vector<MOMAdata> cells =  getData(infile, 
                                            config.time_col,
                                            config.length_col,
                                            config.fp_col,
                                            config.delm,
                                            config.cell_tags,
                                            config.parent_tags);

    /* genealogy built via the parent_id (string) given in data file */
    build_cell_genealogy(cells);
    /* init mean and cov of all root cells */
    init_cells(cells, 5);

    std::cout << "-> Minimizaton" << "\n";
    /* minimization for tree starting from cells[0] */
    minimize_wrapper(&total_likelihood, cells[0], params);
}


int main(int argc, char** argv){
    /* read configuration files */
    CSVconfig config("csv_config.txt");
    std::cout << config;

    Parameter_set params("parameter_bounds.txt");
    std::cout << params;

    /* get input file */
    std::string infile = argv[1];    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }

    /* set and setup (global) output file, */
    std::cout << "Outfile: " << outfile_name(infile, params) << "\n";
    _outfile = outfile_name(infile, params); 
    setup_outfile(_outfile, params);

    /* actual minimization */
    run_minimization(config, params, infile);
    std::cout << "Done." << std::endl;
    return 0;
}
