#include "CSVconfig.h"

#include "likelihood.h"
#include "minimizer_nlopt.h"

#include "tests.h"

#include <iostream> 
#include <iterator> 


void run_minimization(std::vector<MOMAdata> &cells, Parameter_set params, double relative_tol){
    std::cout << "-> Minimizaton" << "\n";
    init_cells(cells, 5);

    /* set and setup (global) output file */
    _outfile = outfile_name_minimization(_infile, params);
    std::cout << "Outfile: " << _outfile << "\n";
    setup_outfile_params(_outfile, params);

    /* minimization for tree starting from cells[0] */
    minimize_wrapper(&total_likelihood, cells[0], params, relative_tol);
}


void run_bound_1dscan(std::vector<MOMAdata> &cells, Parameter_set params){
    std::cout << "-> 1d Scan" << "\n";
    init_cells(cells, 5);
    
    for(int i=0; i<params.all.size(); ++i){
        /* reset params */
        std::vector<double> params_vec = params.get_init();

        if (params.all[i].bound){
            /* 
            * set and setup new (global) output file, for each scan, 
            * filename containing the parameter name that is varied 
            */

            _outfile = outfile_name_scan(_infile, params.all[i].name);
            setup_outfile_params(_outfile, params);
            std::cout << "Outfile: " << _outfile << "\n";

            /* set sampling vector np.arange style*/
            std::vector<double> sampling = arange<double>(params.all[i].lower, 
                                                            params.all[i].upper, 
                                                            params.all[i].step);
            for(int j=0; j<sampling.size(); ++j){
                /* reset mean, cov */
                params_vec[i] = sampling[j];
                total_likelihood(params_vec, cells[0]);
            }
        }
    }
}


void run_prediction(std::vector<MOMAdata> &cells, Parameter_set params){
    std::cout << "-> prediction" << "\n";
    
    _outfile = outfile_name_prediction(_infile);
    std::cout << "Outfile: " << _outfile << "\n";
    setup_outfile_mean(_outfile, params);

    std::vector<double> params_vec = params.get_init();

    /* forward...*/
    init_cells(cells, 5);
    prediction_forward(params_vec, cells[0]);

    /* backward...*/
    init_cells_r(cells, 5);


    /* save */
    std::ofstream file(_outfile, std::ios_base::app);

    for(int i=0; i<cells.size();++i){
        for (int j=0; j<cells[i].mean_forward.size();++j )
            file    << cells[i].mean_forward[j][0] <<','
                    << cells[i].mean_forward[j][1] <<','
                    << cells[i].mean_forward[j][2] <<','
                    << cells[i].mean_forward[j][3] << "\n";  
    }
    file.close();
}


std::map<std::string, std::string> arg_parser(int argc, char** argv){
    std::vector<std::vector<std::string>> keys = {
        {"-h","--help", "\t\thelp message"},
        {"-i", "--infile", "(required) input/data file"},
        {"-b", "--parameter_bounds", "\t(required) file defining the type, step, bounds of the parameters"},
        {"-c", "--csv_config", "\tfile that sets the colums that will be used from the input file"},
        {"-l","--print_level", "\tprint level >=0, default=0 "},
        {"-r","--rel_tol", "\t\trelative tolerance of minimization, default=1e-2"},
        {"-m","--minimize", "\t\trun minimization"},
        {"-s","--scan", "\t\trun 1d parameter scan"},
        {"-p","--predict", "\t\trun prediction"}
        };

    std::map<std::string, int> key_indices; 
    for (int i = 0; i < keys.size(); ++i){
        key_indices.insert(std::pair<std::string, int>(keys[i][0], i)); 
    }

    std::map<std::string, std::string> arguments;
    /* defaults: */
    arguments["csv_config"] = "csv_config.txt";
    arguments["print_level"] = "0";

    arguments["rel_tol"] = "1e-2";

    for(int k=0; k<keys.size(); ++k){
        for(int i=1; i<argc ; ++i){
            if (argv[i] == keys[k][0] || argv[i] == keys[k][1]){
                 if(k==key_indices["-i"]) 
                    arguments["infile"] = argv[i+1];
                else if(k==key_indices["-b"]) 
                    arguments["parameter_bounds"] = argv[i+1];
                else if(k==key_indices["-c"]) 
                    arguments["csv_config"] = argv[i+1];
				else if(k==key_indices["-l"])
                    arguments["print_level"] = argv[i+1];
				else if(k==key_indices["-r"])
                    arguments["rel_tol"] = argv[i+1];
                else if(k==key_indices["-m"])
                    arguments["minimize"] = "1";
                else if(k==key_indices["-s"])
                    arguments["scan"] = "1";
                else if(k==key_indices["-p"])
                    arguments["predict"] = "1";
                else if (k==key_indices["-h"]){
                    arguments["quit"] = "1";
                    std::cout << "Usage: ./gfp_gaussian <infile> [-options]\n";
                    for(int j=0; j<keys.size(); ++j)
                        std::cout << keys[j][0] <<", "<< keys[j][1] << keys[j][2] <<"\n";
                }
            }
        }
    }
    /* Check is required filenames are parsed and files exist */
    if (!arguments.count("infile")){
        std::cout << "Required infile flag not set!\n";
        arguments["quit"] = "1";
    }
    else if(! std::__fs::filesystem::exists(arguments["infile"])){
        std::cout << "Infile " << _infile << " not found (use '-h' for help)!" << std::endl;
        arguments["quit"] = "1";
    }

    if (!arguments.count("parameter_bounds")){
        std::cout << "Required parameter_bounds flag not set!\n";
        arguments["quit"] = "1";
    }
    else if(! std::__fs::filesystem::exists(arguments["parameter_bounds"])){   
        std::cout << "Paramters bound file " << arguments["parameter_bounds"] << " not found (use '-h' for help)!" << std::endl;
        arguments["quit"] = "1";
    }

    /* Check if csv file (if parsed) exists, to avoid confusion */
    if(arguments.count("csv_config") && !std::__fs::filesystem::exists(arguments["csv_config"])){   
        std::cout << "csv_config flag set, but csv configuration file " << arguments["csv_config"] << " not found!" << std::endl;
        arguments["quit"] = "1";
    }
    return arguments;
}


int main(int argc, char** argv){
    /* read configuration files */
    std::map<std::string, std::string> arguments = arg_parser(argc, argv);
    _print_level = std::stoi(arguments["print_level"]);

    if (arguments.count("quit")){
        std::cout << "Quit\n";
        return 0;    
    }

    /* get parameter and csv config file */
    Parameter_set params(arguments["parameter_bounds"]);
    std::cout << params << "\n";

    CSVconfig config(arguments["csv_config"]);
    std::cout << config << "\n";

    /* Read data from input file */
    _infile = arguments["infile"];    
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

    /* bound_1dscan, minimization... */
    if (arguments.count("minimize"))
        run_minimization(cells, params, std::stod(arguments["rel_tol"]));

    else if (arguments.count("scan"))
        run_bound_1dscan(cells, params);

    else if (arguments.count("predict"))
        run_prediction(cells, params);

    std::cout << "Done." << std::endl;
    return 0;
}


