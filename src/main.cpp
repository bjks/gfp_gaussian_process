#include "CSVconfig.h"
#include "likelihood.h"
#include "minimizer_nlopt.h"
#include "tests.h"

#include <iostream> 
#include <iterator> 


void run_minimization(std::vector<MOMAdata> &cells, Parameter_set params, double relative_tol){
    std::cout << "-> Minimizaton" << "\n";

    /* set and setup (global) output file */
    _outfile = outfile_name_minimization(_infile, params);
    std::cout << "Outfile: " << _outfile << "\n";
    setup_outfile(_outfile, params);

    /* minimization for tree starting from cells[0] */
    minimize_wrapper(&total_likelihood, cells[0], params, relative_tol);
}


void bound_1dscan(std::vector<MOMAdata> &cells, Parameter_set params){
    std::cout << "-> 1d Scan" << "\n";
    
    for(int i=0; i<params.all.size(); ++i){
        /* reset params */
        std::vector<double> params_vec = params.get_init();

        if (params.all[i].bound){
            /* 
            * set and setup new (global) output file, for each scan, 
            * filename containing the parameter name that is varied 
            */

            _outfile = outfile_name_scan(_infile, params.all[i].name);
            setup_outfile(_outfile, params);
            std::cout << "Outfile: " << _outfile << "\n";

            /* set sampling vector np.arange style*/
            std::vector<double> sampling = arange<double>(params.all[i].lower, 
                                                            params.all[i].upper, 
                                                            params.all[i].step);
            for(int j=0; j<sampling.size(); ++j){
                /* reset mean, cov */
                params_vec[i] = sampling[j];
                total_likelihood(params_vec, &cells[0]);
            }
        }
    }
}



std::map<std::string, std::string> arg_parser(int argc, char** argv){
    std::vector<std::vector<std::string>> keys = {{"-p", "--parameter_config", "\tfile defining the type, step, bounds of the parameters"},
                                                {"-c", "--csv_config", "\tfile that sets the colums that will be used from the input file"},
                                                {"-m","--mode", "\t\tmode keyword can start with 'm'->minimization or 's'->scan"},
                                                {"-l","--print_level", "\tprint level >=0 "},
                                                {"-r","--rel_tol", "\t\trelative tolerance of minimization"},
                                                {"-h","--help", "\t\thelp message"}};

    std::map<std::string, std::string> arguments;
    /* defaults: */
    arguments["parameter_config"] = "parameter_bounds.txt";
    arguments["csv_config"] = "csv_config.txt";
    arguments["mode"] = "minimize";
    arguments["rel_tol"] = "1e-2";

    for(int k=0; k<keys.size(); ++k){
        for(int i=1; i<argc ; ++i){
            if (argv[i] == keys[k][0] || argv[i] == keys[k][1]){
                if(k==0) 
                    arguments["parameter_config"] = argv[i+1];
                else if(k==1) 
                    arguments["csv_config"] = argv[i+1];
				else if(k==2) 
                    arguments["mode"] = argv[i+1];
				else if(k==3)
                    _print_level = std::stoi(argv[i+1]);
				else if(k==4)
                    arguments["rel_tol"] = argv[i+1];
                else if (k==5){
                    arguments["mode"] = "quit";
                    std::cout << "Usage: ./gfp_gaussian <infile> [-options]\n";
                    for(int j=0; j<keys.size(); ++j)
                        std::cout << keys[j][0] <<", "<< keys[j][1] << keys[j][2] <<"\n";
                }
            }
        }
    }
    return arguments;
}


int main(int argc, char** argv){
    /* read configuration files */
    std::map<std::string, std::string> arguments = arg_parser(argc, argv);

    if (arguments["mode"] == "quit")
        return 0;    

    /* get input file */
    _infile = argv[1];    
    if(! std::__fs::filesystem::exists(_infile)){
        std::cout << "File " << _infile << " not found (use '-h' for help)! \nQuit" << std::endl;
        return 0;
    }

    /* get parameter file */
    if(! std::__fs::filesystem::exists(arguments["parameter_config"])){   
        std::cout << "File " << arguments["parameter_config"] << " not found (use '-h' for help)! \nQuit" << std::endl;
        return 0;
    }
    Parameter_set params(arguments["parameter_config"]);
    std::cout << params;

    /* get csv config file */
    if(! std::__fs::filesystem::exists(arguments["csv_config"])){   
        std::cout << "File " << arguments["parameter_config"] << " not found! \nUse defaults" << std::endl;
    }
    CSVconfig config(arguments["csv_config"]);
    std::cout << config;


    /* Read data from input file */
    std::cout << "-> Reading" << "\n";
    std::vector<MOMAdata> cells =  getData(_infile, 
                                            config.time_col,
                                            config.length_col,
                                            config.fp_col,
                                            config.delm,
                                            config.cell_tags,
                                            config.parent_tags);

    /* genealogy built via the parent_id (string) given in data file */
    build_cell_genealogy(cells);
    init_cells(cells, 5);

    /* bound_1dscan, run_minimization... */
    if (arguments["mode"].at(0) == 'm')
        run_minimization(cells, params, std::stod(arguments["rel_tol"]));

    else if (arguments["mode"].at(0) == 's')
        bound_1dscan(cells, params);

    std::cout << "Done." << std::endl;
    return 0;
}


