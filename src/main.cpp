#include "likelihood.h"
#include "minimizer_nlopt.h"

#include "tests.h"
#include <filesystem>
#include <iostream> 
#include <iterator> 
#include <iomanip> 


int run_minimization(std::vector<MOMAdata> &cells, 
                    Parameter_set &params, 
                    std::map<std::string, std::string> arguments){

    std::cout << "-> Minimizaton" << "\n";

    /* set and setup (global) output file */
    _outfile_ll = outfile_name_minimization_process(arguments, params);
    setup_outfile_likelihood(_outfile_ll, params);
    std::cout << "Outfile: " << _outfile_ll << "\n";

    /* minimization for tree starting from cells[0] */
    bool found_min = false;
    std::string min_algo = "LN_COBYLA";
    double ll_max;
    if (arguments["search_space"] == "log"){
        ll_max = minimize_wrapper_log_params(&total_likelihood_log_params, cells, 
                                        params, std::stod(arguments["tolerance"]), 
                                        found_min, min_algo);
    }
    else{
        ll_max = minimize_wrapper(&total_likelihood, cells, 
                                        params, std::stod(arguments["tolerance"]), 
                                        found_min, min_algo);
    }

    if (!found_min){
        return -1;
    }
    /* estimate errors of params via hessian */
    std::cout << "-> Error estimation" << "\n";
    std::string outfile_estim = outfile_name_minimization_final(arguments, params);
    std::cout << "Outfile: " << outfile_estim << "\n";

    params.to_csv(outfile_estim);
    save_error_bars(outfile_estim, params, cells);
    save_final_likelihood(outfile_estim, cells, ll_max, min_algo, 
                        std::stod(arguments["tolerance"]), 
                        arguments["search_space"]);

    return 0;
}


void run_bound_1dscan(std::vector<MOMAdata> &cells, Parameter_set params,
                      std::map<std::string, std::string> arguments){
    std::cout << "-> 1d Scan" << "\n";
    _save_ll = true;
    for(size_t i=0; i<params.all.size(); ++i){
        if (params.all[i].bound){
            /* 
            reset params, if paramters have been minimized the final value will be taken
            otherwise the init values are taken 
            */
            std::vector<double> params_vec = params.get_final();
             
            /* 
            set and setup new (global) output file, for each scan, 
            filename containing the parameter name that is varied
            */
            _outfile_ll = outfile_name_scan(arguments, params.all[i].name);
            setup_outfile_likelihood(_outfile_ll, params);
            std::cout << "Outfile: " << _outfile_ll << "\n";

            /* set sampling vector np.arange style*/
            std::vector<double> sampling = arange<double>(params.all[i].lower, 
                                                            params.all[i].upper, 
                                                            params.all[i].step);
            for(size_t j=0; j<sampling.size(); ++j){
                params_vec[i] = sampling[j];
                total_likelihood(params_vec, cells);
            }
        }
    }
    _save_ll = false;
}


void run_prediction(std::vector<MOMAdata> &cells, Parameter_set params, 
                    std::map<std::string, std::string> arguments, const CSVconfig &config){
    std::cout << "-> prediction" << "\n";
    
    std::string outfile = outfile_name_prediction(arguments);
    std::string outfile_b = outfile_name_prediction(arguments, "_backward");
    std::string outfile_f = outfile_name_prediction(arguments, "_forward");


    std::cout << "Outfile: " << outfile << "\n";
    std::cout << "Outfile forward: " << outfile_f << "\n";
    std::cout << "Outfile backward: " << outfile_b << "\n";


    std::vector<double> params_vec = params.get_final();

    /* forward...*/
    prediction_forward(params_vec, cells);

    /* backward...*/
    prediction_backward(params_vec, cells);

    /* combine the two */
    combine_predictions(cells);

    /* save */
    write_pretictions_to_file(cells, outfile_b, params, config, "b");
    write_pretictions_to_file(cells, outfile_f, params, config, "f");

    write_pretictions_to_file(cells, outfile, params, config);
}


std::map<std::string, std::string> arg_parser(int argc, char** argv){
    std::vector<std::vector<std::string>> keys = {
        {"-h","--help", "this help message"},
        {"-i", "--infile", "(required) input data file"},
        {"-b", "--parameter_bounds", "(required) file setting the type, step, bounds of the parameters"},
        {"-c", "--csv_config", "file that sets the colums that will be used from the input file"},
        {"-l","--print_level", "print level >=0, default=0"},
        {"-o","--outdir", "specify output direction and do not use default"},
        {"-t","--tolerance", "absolute tolerance of maximization between optimization steps, default=1e-1"},
        {"-space","--search_space", "search parameter space in 'log' or 'linear' space, default: 'linear'"},
        {"-stat","--stationary", "indicates that the cell are not growing (much)"},
        {"-m","--maximize", "run maximization"},
        {"-s","--scan", "run 1d parameter scan"},
        {"-p","--predict", "run prediction"}
        };

    std::map<std::string, int> key_indices; 
    for (size_t i = 0; i < keys.size(); ++i){
        key_indices.insert(std::pair<std::string, int>(keys[i][0], i)); 
    }

    std::map<std::string, std::string> arguments;
    /* defaults: */
    arguments["print_level"] = "0";
    arguments["tolerance"] = "1e-1";
    arguments["search_space"] = "linear";


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
				else if(k==key_indices["-o"])
                    arguments["outdir"] = argv[i+1];
				else if(k==key_indices["-t"])
                    arguments["tolerance"] = argv[i+1];
                else if(k==key_indices["-space"])
                    arguments["search_space"] = argv[i+1];
                else if(k==key_indices["-stat"])
                    arguments["stationary"] = "1";
                else if(k==key_indices["-m"])
                    arguments["minimize"] = "1";
                else if(k==key_indices["-s"])
                    arguments["scan"] = "1";
                else if(k==key_indices["-p"])
                    arguments["predict"] = "1";
                else if (k==key_indices["-h"]){
                    arguments["quit"] = "1";
                    std::cout << "Usage: ./gfp_gaussian [-options]\n";
                    for(size_t j=0; j<keys.size(); ++j)
                        std::cout << pad_str(keys[j][0] + ", "+ keys[j][1], 27) << keys[j][2] <<"\n";
                }
            }
        }
    }
    if (arguments.count("quit")){
        return arguments;
    }

    /* Check if meaningfull search space argument */
    if (arguments["search_space"] != "log" && arguments["search_space"] != "linear"){
        std::cerr << "search_space must be either 'log' or 'linear', not " << arguments["search_space"];
        arguments["quit"] = "1";
    }

    /* Check is required filenames are parsed and files exist */
    if (!arguments.count("infile")){
        std::cerr << "Required infile flag not set!\n";
        arguments["quit"] = "1";
    }
    else if(! std::filesystem::exists(arguments["infile"])){
        std::cerr << "Infile " << arguments["infile"] << " not found (use '-h' for help)!" << std::endl;
        arguments["quit"] = "1";
    }

    if (!arguments.count("parameter_bounds")){
        std::cerr << "Required parameter_bounds flag not set!\n";
        arguments["quit"] = "1";
    }
    else if(! std::filesystem::exists(arguments["parameter_bounds"])){   
        std::cerr << "Paramters bound file " << arguments["parameter_bounds"] << " not found (use '-h' for help)!" << std::endl;
        arguments["quit"] = "1";
    }

    /* Check if csv file (if parsed) exists, to avoid confusion */
    if(arguments.count("csv_config") && !std::filesystem::exists(arguments["csv_config"])){   
        std::cerr << "csv_config flag set, but csv configuration file " << arguments["csv_config"] << " not found!" << std::endl;
        arguments["quit"] = "1";
    }
    return arguments;
}


int main(int argc, char** argv){
    /* process command line arguments */
    std::map<std::string, std::string> arguments = arg_parser(argc, argv);
    _print_level = std::stoi(arguments["print_level"]);

    if (arguments.count("quit")){
        std::cout << "Quit\n";
        return 0;    
    }

    /* get parameter and csv config file */
    Parameter_set params(arguments["parameter_bounds"]);
    if (!params.is_complete()){
        std::cout << "Quit\n";
        return -1;    
    }
    std::cout << params << "\n";

    CSVconfig config(arguments["csv_config"]);
    std::cout << config << "\n";

    /* Read data from input file */
    std::cout << "-> Reading" << "\n";
    std::vector<MOMAdata> cells =  get_data(arguments["infile"], config);
    if (!cells.size()){
        std::cout << "Quit\n";
        return -1;    
    }
    /* genealogy built via the parent_id (string) given in data file */
    bool valid_genealogy = build_cell_genealogy(cells);
    if (!valid_genealogy){
        return -1;
    }

    /* inititialize mean and cov for forward and backward direction */
    init_cells(cells, arguments.count("stationary"));
    init_cells_r(cells, arguments.count("stationary"));


    /* run bound_1dscan, minimization and/or prediction... */
    if (arguments.count("minimize")){
        int min_message = run_minimization(cells, params, arguments);
        if (min_message  == -1){
            return -1;
        }
    }
    if (arguments.count("scan"))
        run_bound_1dscan(cells, params, arguments);

    if (arguments.count("predict"))
        run_prediction(cells, params, arguments, config);

    std::cout << "Done." << std::endl;
    return 0;
}


