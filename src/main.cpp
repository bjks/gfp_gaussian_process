#include "CSVconfig.h"
#include "likelihood.h"
#include "minimizer_nlopt.h"

#include <iostream> 
#include <iterator> 



int main(int argc, char** argv){
    // Eigen::VectorXd v(3);
    // Eigen::VectorXd w(3);

    // v << -1, -2,  -3;
    // w << 10,  20, 30;

    // Eigen::VectorXd x(2);
    // x << 2, 4;
    // std::cout << x  << std::endl;

    // Eigen::MatrixXd m(2, 3);

    // m << v.transpose(), w.transpose();


    // std::cout << m  << std::endl;
    // std::cout << v  << std::endl;
    // std::cout << rowwise_add(m, x)  << std::endl;

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
    return 0;


    std::cout << "-> Minimizaton" << "\n";
    // minimization for tree starting from cells[0]
    minimize_wrapper(&total_likelihood, cells[0], params);


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
