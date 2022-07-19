#include <fstream>
#include <string>

#include <filesystem>
#include <iostream> 
#include <iterator> 
#include <iomanip> 

std::ofstream _file_log;

#include "likelihood.h"
#include "minimizer_nlopt.h"
// #include "tests.h"

// #include <boost/iostreams/filtering_streambuf.hpp>
// #include <boost/iostreams/copy.hpp>
// #include <boost/iostreams/filter/gzip.hpp>

/* ======================================== */
/* Routines for the different running modes */
/* ======================================== */

void run_minimization(std::vector<MOMAdata> &cells, 
                    Parameter_set &params, 
                    std::map<std::string, 
                    std::string> arguments, 
                    int segment){

    _file_log << "-> Minimizaton" << "\n";

    /* set and setup (global) output file */
    std::string outfile_ll = outfile_name_minimization_process(arguments, params, segment);
    _file_iteration = std::ofstream(outfile_ll, std::ios_base::app);

    setup_outfile_likelihood(outfile_ll, params);
    _file_log << "Outfile: " << outfile_ll << "\n";

    /* minimization for tree starting from cells[0] */
    std::string min_algo = "LN_NELDERMEAD";
    double ll_max;
    if (arguments["search_space"] == "log"){
        ll_max = minimize_wrapper_log_params(&total_likelihood_log_params, cells, 
                                        params, std::stod(arguments["tolerance_maximization"]), 
                                        min_algo);
    }
    else{
        ll_max = minimize_wrapper(&total_likelihood, cells, 
                                        params, std::stod(arguments["tolerance_maximization"]), 
                                        min_algo);
    }

    _file_iteration.close();

    /* estimate errors of params via hessian */
    _file_log << "-> Error estimation" << "\n";
    std::string outfile_estim = outfile_name_minimization_final(arguments, params, segment);
    _file_log << "Outfile: " << outfile_estim << "\n";

    params.to_csv(outfile_estim);
    save_error_bars(outfile_estim, params, cells);
    save_final_likelihood(outfile_estim, 
                          cells, 
                          ll_max, 
                          min_algo, 
                          std::stod(arguments["tolerance_maximization"]), 
                          arguments["search_space"], 
                          arguments["noise_model"],
                          arguments["cell_division_model"],
                          "0.4.1");

    std::string outfile_params = outfile_name_parameter_file(arguments, params, segment);
    create_parameter_file(outfile_params, params);

}


void run_bound_1dscan(std::vector<MOMAdata> &cells, 
                    Parameter_set params,
                    std::map<std::string, std::string> arguments, 
                    int segment){
    _file_log << "-> 1d Scan" << "\n";
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
            std::string outfile_scan = outfile_name_scan(arguments, params.all[i].name, segment);
            _file_iteration = std::ofstream(outfile_scan, std::ios_base::app);

            setup_outfile_likelihood(outfile_scan, params);
            _file_log << "Outfile: " << outfile_scan << "\n";

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


void run_prediction_segments(std::vector<MOMAdata> &cells, 
                            std::vector<Parameter_set> params_list, 
                            std::map<std::string, std::string> arguments,
                            const CSVconfig &config){
    _file_log << "-> prediction" << "\n";
    
    std::string outfile   = outfile_name_prediction(arguments, params_list);
    std::string outfile_b = outfile_name_prediction(arguments, params_list, "_backward");
    std::string outfile_f = outfile_name_prediction(arguments, params_list, "_forward");

    std::vector<std::vector<double>> params_vecs;
    for (size_t i=0; i<params_list.size(); ++i){
        params_vecs.push_back(params_list[i].get_final());
    }

    // forward...
    _file_log << "Outfile forward: " << outfile_f << "\n";
    prediction_forward(params_vecs, cells);

    // backward...
    _file_log << "Outfile backward: " << outfile_b << "\n";
    prediction_backward(params_vecs, cells);

    // combine the two 
    _file_log << "Outfile: " << outfile << "\n";
    combine_predictions(cells);

    /* save */
    write_predictions_to_file(cells, outfile_f, params_list, "f");
    write_predictions_to_file(cells, outfile_b, params_list, "b");

    write_predictions_to_file(cells, outfile, params_list);
}


void run_joint_distribution(std::vector<MOMAdata> &cells, 
                    std::vector<Parameter_set> params_list, 
                    std::map<std::string, std::string> arguments, 
                    const CSVconfig &config){
    _file_log << "-> joint posteriors" << "\n";

    std::vector<std::vector<double>> params_vecs;
    for (size_t i=0; i<params_list.size(); ++i){
        params_vecs.push_back(params_list[i].get_final());
    }

    std::string outfile_joints = outfile_name_joints(arguments, params_list);

    /* if I get Boost to work on the cluster that is an option */
    // std::ofstream file(outfile_joints + ".gz", std::ios_base::out | std::ios_base::binary);
    // boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
    // outbuf.push(boost::iostreams::gzip_compressor());
    // outbuf.push(file);
    // //Convert streambuf to ostream
    // std::ostream out(&outbuf);

    // setup_outfile_joints(out, params_list);

    // // calculate all possible joints
    // collect_joint_distributions(params_vecs, cells, out);

    // boost::iostreams::close(outbuf);
    // file.close();

    std::ofstream file(outfile_joints);
    setup_outfile_joints(file, params_list);

    // calculate all possible joints
    collect_joint_distributions(params_vecs, cells, file, std::stod(arguments["rel_tolerance_joints"]));
    file.close();

}



/* =============================================================================== */
std::map<std::string, std::string> arg_parser(int argc, char** argv){
    std::vector<std::vector<std::string>> keys = 
        {
        {"-h",      "--help",                   "this help message"},
        {"-i",      "--infile",                 "(required) input data file"},
        {"-b",      "--parameter_bounds",       "(required) file(s) setting the type, step, bounds of the parameters"},
        {"-c",      "--csv_config",             "file that sets the colums that will be used from the input file"},
        {"-l",      "--print_level",            "print level {0,1,2}, default: 0"},
        {"-o",      "--outdir",                 "specify output direction and do not use default"},
        {"-t",      "--tolerance_maximization",  "absolute tolerance of maximization between optimization steps, default: 1e-10"},
        {"-r",      "--rel_tolerance_joints",    "relative tolerance of joint calculation: default 1e-10"},
        {"-space",  "--search_space",           "search parameter space in {'log'|'linear'} space, default: 'log'"},
        {"-noise",  "--noise_model",            "measurement noise of fp content {'scaled'|'const'} default: 'scaled'"},
        {"-div",    "--cell_division_model",          "cell divison model {'binomial'|'gauss'} default: 'binomial'"},
        {"-m",      "--maximize",               "run maximization"},
        {"-s",      "--scan",                   "run 1d parameter scan"},
        {"-p",      "--predict",                "run prediction"},
        {"-j",      "--joints",                 "run calculation of joint probabilities"}
        };

    std::map<std::string, int> key_indices; 
    for (size_t i = 0; i < keys.size(); ++i){
        key_indices.insert(std::pair<std::string, int>(keys[i][0], i)); 
    }

    std::map<std::string, std::string> arguments;
    /* defaults: */
    arguments["print_level"] = "0";
    arguments["tolerance_maximization"] = "1e-10";
    arguments["rel_tolerance_joints"] = "1e-10";
    arguments["search_space"] = "log";
    arguments["noise_model"] = "scaled";
    arguments["cell_division_model"] = "binomial";


    for(int k=0; k<keys.size(); ++k){
        for(int i=1; i<argc ; ++i){
            if (argv[i] == keys[k][0] || argv[i] == keys[k][1]){
                 if(k==key_indices["-i"]) 
                    arguments["infile"] = argv[i+1];
                else if(k==key_indices["-b"]){
                    for(int j=i+1; j<argc ; ++j){
                        std::string argj = argv[j];
                        if (argj.rfind("-", 0) == 0 ){
                            break;
                        }
                        arguments["parameter_bounds"] += argj + " ";
                    }
                    arguments["parameter_bounds"] = trim(arguments["parameter_bounds"]);
                }
                else if(k==key_indices["-c"]) 
                    arguments["csv_config"] = argv[i+1];
				else if(k==key_indices["-l"])
                    arguments["print_level"] = argv[i+1];
				else if(k==key_indices["-o"])
                    arguments["outdir"] = argv[i+1];
				else if(k==key_indices["-t"])
                    arguments["tolerance_maximization"] = argv[i+1];
                else if(k==key_indices["-r"])
                    arguments["rel_tolerance_joints"] = argv[i+1];
                else if(k==key_indices["-space"])
                    arguments["search_space"] = argv[i+1];
                else if(k==key_indices["-noise"])
                    arguments["noise_model"] = argv[i+1];
                else if(k==key_indices["-div"])
                    arguments["cell_division_model"] = argv[i+1];
                else if(k==key_indices["-space"])
                    arguments["search_space"] = argv[i+1];
                else if(k==key_indices["-m"])
                    arguments["minimize"] = "1";
                else if(k==key_indices["-s"])
                    arguments["scan"] = "1";
                else if(k==key_indices["-p"])
                    arguments["predict"] = "1";
                else if(k==key_indices["-j"]){
                    arguments["joints"] = "1";
                    arguments["predict"] = "1"; // needs to be run prior to the auto covariance calculation
                }
                else if (k==key_indices["-h"]){
                    arguments["help"] = "1";
                    _file_log << "Usage: ./gfp_gaussian [-options]\n";
                    for(size_t j=0; j<keys.size(); ++j)
                        _file_log << pad_str(keys[j][0] + ", "+ keys[j][1], 35) << keys[j][2] <<"\n";
                }
            }
        }
    }
    if (arguments.count("help")){
        return arguments;
    }

    /* Check if meaningfull search space argument */
    if (arguments["search_space"] != "log" && arguments["search_space"] != "linear"){
        _file_log << "(arg_parser) ERROR: search_space must be either 'log' or 'linear', not " << arguments["search_space"];
        throw std::invalid_argument("Invalide argument");
    }

    if (arguments["noise_model"] != "const" && arguments["noise_model"] != "scaled"){
        _file_log << "(arg_parser) ERROR: noise_model must be either 'const' or 'scaled', not " << arguments["noise_model"];
        throw std::invalid_argument("Invalide argument");
    }
    if (arguments["cell_division_model"] != "gauss" && arguments["cell_division_model"] != "binomial"){
        _file_log << "(arg_parser) ERROR: cell_division_model must be either 'gauss' or 'binomial', not " << arguments["cell_division_model"];
        throw std::invalid_argument("Invalide argument");
    }

    /* Check is required filenames are parsed and files exist */
    if (!arguments.count("infile")){
        _file_log << "(arg_parser) ERROR: Required infile flag not set!\n";
        throw std::invalid_argument("Invalide argument");
    }
    else if(! std::filesystem::exists(arguments["infile"])){
        _file_log << "(arg_parser) ERROR: Infile " << arguments["infile"] << " not found (use '-h' for help)!" << std::endl;
        throw std::invalid_argument("Invalide argument");
    }

    if (!arguments.count("parameter_bounds")){
        _file_log << "(arg_parser) ERROR: Required parameter_bounds flag not set!\n";
        throw std::invalid_argument("Invalide argument");
    }


    std::vector<std::string> param_files = split_string_at(arguments["parameter_bounds"], " ");
    for (size_t i=0; i<param_files.size(); ++i){
        if(!std::filesystem::exists(param_files[i])){   
            _file_log << "(arg_parser) ERROR: Paramters bound file '" << param_files[i] << "' not found (use '-h' for help)!" << std::endl;
            throw std::invalid_argument("Invalide argument");
        }
    }

    /* Check if csv file (if parsed) exists, to avoid confusion */
    if(arguments.count("csv_config") && !std::filesystem::exists(arguments["csv_config"])){   
        _file_log << "(arg_parser) ERROR: csv_config flag set, but csv configuration file " << arguments["csv_config"] << " not found!" << std::endl;
        throw std::invalid_argument("Invalide argument");
    }
    return arguments;
}

std::string outfile_name_log(std::map<std::string, std::string> arguments, 
                                    std::vector<Parameter_set> params_list, std::string suffix=""){
    /* Filename for a log file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]);
    for (size_t i=0; i<params_list.size(); ++i){
        outfile += outfile_param_code(params_list[i]);
    }
    return outfile + suffix + ".log";
}

/* =============================================================================== */
/*                                  MAIN                                           */
/* =============================================================================== */

int main(int argc, char** argv){

    // test_mean_cov_model();
    // return 0;
    std::string outfile_log;
    std::string outfile_log_success;
    std::string outfile_log_error;

    std::cout << "Running... \n";

    try{
        /* process command line arguments */
        std::map<std::string, std::string> arguments = arg_parser(argc, argv);
        _print_level = std::stoi(arguments["print_level"]);

        if (arguments.count("help")){
            return EXIT_SUCCESS;    
        }

        /* Read parameters as a vector of Parameter_set instances */
        std::vector<std::string> param_files = split_string_at(arguments["parameter_bounds"], " ");
        std::vector<Parameter_set> params_list;
        
        outfile_log  = outfile_name_log(arguments, params_list);
        outfile_log_success  = outfile_name_log(arguments, params_list, "_success");
        outfile_log_error  = outfile_name_log(arguments, params_list, "_error");
        
        _file_log = std::ofstream(outfile_log, std::ios_base::app);
        std::cout << "Temporary log file '" << outfile_log <<  "' created\n";

        for (size_t i=0; i<param_files.size(); ++i){
            Parameter_set params(param_files[i]);
            params.check_if_complete();
            _file_log << params << "\n";

            params_list.push_back(params);
        }

        /* Read csv config file */
        CSVconfig config(arguments["csv_config"]);
        _file_log << config << "\n";

        /* Read data from input file */
        _file_log << "-> Reading" << "\n";
        std::vector<MOMAdata> cells =  read_data(arguments["infile"], 
                                                config, 
                                                arguments["noise_model"],
                                                arguments["cell_division_model"]);

        /* Count segments and check if we have enough parameter files */
        std::vector<int> segment_indices = get_segment_indices(cells);
        if (segment_indices.size() != param_files.size()){
            _file_log   << "(main) ERROR: There are " << segment_indices.size() 
                        << " segments, but " << param_files.size() << " parameter files!\n";
            throw std::invalid_argument("Invalide argument");
        }


        /* ============================================================== */
        /* run bound_1dscan, minimization and/or prediction... */

        if (arguments.count("minimize")){
            for(size_t i=0; i<segment_indices.size(); ++i){
                if (params_list[i].has_nonfixed()){
                    std::vector<MOMAdata> cells_in_segment = get_segment(cells, segment_indices[i]);

                    /* genealogy built via the parent_id (string) given in data file */
                    build_cell_genealogy(cells_in_segment);

                    /* inititialize mean and cov for forward and backward direction */
                    init_cells(cells_in_segment);

                    /* Run actual minimization */
                    run_minimization(cells_in_segment, params_list[i], arguments, get_segment_file_number(segment_indices, i));
                }
            }
        }

        if (arguments.count("scan")){
            for(size_t i=0; i<segment_indices.size(); ++i){
                std::vector<MOMAdata> cells_in_segment = get_segment(cells, segment_indices[i]);

                /* genealogy built via the parent_id (string) given in data file */
                build_cell_genealogy(cells_in_segment);

                /* inititialize mean and cov for forward and backward direction */
                init_cells(cells_in_segment);

                /* Run scan*/
                run_bound_1dscan(cells_in_segment, params_list[i], arguments, get_segment_file_number(segment_indices, i));
            }
        }

        if (arguments.count("predict")){
            build_cell_genealogy(cells);
            init_cells(cells);
            run_prediction_segments(cells, params_list, arguments, config);
        }
        
        if (arguments.count("joints")){
            /* Run joints, note that prediction was already run as this stage! */
            run_joint_distribution(cells, params_list, arguments, config);
        }

        _file_log << "Done." << std::endl;

        std::cout << "Done. Log file: " << outfile_log_success << std::endl;

        std::rename(outfile_log.c_str(), outfile_log_success.c_str());
        _file_log.close();

        return EXIT_SUCCESS;
    }
    catch (std::exception &e) {
        std::rename(outfile_log.c_str(), outfile_log_error.c_str());
        _file_log << "Quit because of an error: " << e.what() <<"\n";
        _file_log.close();

        std::cout << "Quit because of an error: " << e.what() << "\n";
        std::cout << "Error log file: " << outfile_log_error << std::endl;

        return EXIT_FAILURE;
    }
}
