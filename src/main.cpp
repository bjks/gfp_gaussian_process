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

    /* genealogy */
    build_cell_genealogy(cells);

    init_cells(cells, 5);
    std::cout << cells[0] << "\n";

    std::cout << "-> Minimizaton" << "\n";
    /* minimization for tree starting from cells[0] */
    minimize_wrapper(&total_likelihood, cells[0], params);
}

void run_likelihood(CSVconfig config, Parameter_set params, std::string infile){

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

    Eigen::VectorXd mean(4);
    mean << 6.93147181e-01,
            6.02556189e+03,
            1.03989065e-02,
            1.02454986e+01;

    Eigen::MatrixXd cov(4,4);
    cov <<  1.22158419e-04, -3.38642002e-01,  3.42444314e-06, -4.90827026e-04,
            -3.38642002e-01,  1.25286734e+05, -3.88680250e-01,  1.42591667e+02,
            3.42444314e-06, -3.88680250e-01,  4.47368172e-06,  5.05127089e-05,
            -4.90827026e-04,  1.42591667e+02,  5.05127089e-05,  2.38674307e+00;

    init_cells(cells, mean, cov);

    double tl = 0;
    pvector(params.get_init());
    sc_likelihood(params.get_init(), cells[0], tl);
    std::cout << "tl: " << tl << "\n";
}

int main(int argc, char** argv){
  
    CSVconfig config("csv_config.txt");
    std::cout << config;

    Parameter_set params("parameter_bounds.txt");
    std::cout << params;

    std::string infile = argv[1];    
    if(! std::__fs::filesystem::exists(infile)){
        std::cout << "File " << infile << " not found! \nQuit" << std::endl;
        return 0;
    }

    run_minimization(config, params, infile);
    std::cout << "Done." << std::endl;
    return 0;
}
