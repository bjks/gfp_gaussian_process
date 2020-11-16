#include "CSVconfig.h"
#include "likelihood.h"
#include "minimizer_nlopt.h"


#include <iostream> 
#include <iterator> 



int main(int argc, char** argv){
    // double x = zerotauint(0.1, 0.1, 0.1, 180, 0);
    // double y = onetauint(0.1, 0.1, 0.1, 20, 0);
    // double z = twotauint(0.1, 0.1, 0.1, 20, 0);
    // double a = treetauint(0.1, 0.1, 0.1, 1e10, 0);



    // std::cout << x << "\n";
    // std::cout << y << "\n";
    // std::cout << z << "\n";
    // std::cout << a << "\n";

    // return 0;

    CSVconfig config("csv_config.txt");
    Parameter_set params("parameter_bounds.txt");

    std::cout << params << "\n";

    // Read data
    std::cout << "-> Reading" << "\n";
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

    // genealogy
    build_cell_genealogy(cells);
    // print_cells(cells);

    init_cells(cells);
    
    MOMAdata cell = cells[0];
    std::cout << cell.cov << "\n" << cell.mean << "\n";
    
    std::vector<double> params_vec = params.get_values();
    double total_likelihood = 0;
    sc_likelihood(params_vec,cell, total_likelihood);
    std::cout << cell.cov << "\n" << cell.mean << "\n";

    return 0 ;
    std::cout << "-> Minimizaton" << "\n";
    // minimization for tree starting from cells[0]
    // minimize_wrapper(&total_likelihood, cells[0], params);


    // // DEMO on how this works for all trees
    // // get the "trees" starting from all root cells
    // std::vector<MOMAdata *> root_cells = get_roots(cells);
    // std::vector<double> dummy_params;
    // for(long j=0; j<root_cells.size(); ++j){
    //     print_generation_tree(dummy_params, *root_cells[j]);
    // }
    
    std::cout << "Done." << std::endl;
    return 0;
}
