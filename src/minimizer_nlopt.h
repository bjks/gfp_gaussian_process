#include <iostream>
#include <vector>
#include <iomanip> 
#include <nlopt.hpp>
#include <math.h>

#include "Parameters.h"


double minimize_wrapper(double (*target_func)(const std::vector<double> &x, std::vector<double> &grad, void *p),
                        std::vector<MOMAdata> &cells,
                        Parameter_set &params, 
                        double tolerance, 
                        std::string opt_name){

    /* 
    * Wraps the minimization step using the nlopt library. 
    * The target function is the log_likelihood calculation returning -log_likelihood.
    * Returns the maximum (!) of the log_likelihood
    * Stopping caused by rounding issues during minimization are caught and the last paramter state is taken.
    */                        

    // set parameter space 
    std::vector<double> lower_bounds(params.all.size());
    std::vector<double> upper_bounds(params.all.size());
    std::vector<double> steps(params.all.size());

    std::vector<double> parameter_state(params.all.size());

    for (size_t i=0; i<params.all.size(); ++i){
        if (params.all[i].minimized){
            parameter_state[i] = params.all[i].final;
        }
        else{
            parameter_state[i] = params.all[i].init;
        }
        if (params.all[i].fixed){
            steps[i] = 1.;                          // will not be used anyway, needs to be non-zero
            lower_bounds[i] = params.all[i].init;
            upper_bounds[i] = params.all[i].init;
        } else{
            steps[i] = params.all[i].step;
            lower_bounds[i] = params.all[i].lower;
            upper_bounds[i] = params.all[i].upper;
        }
    }

    // set up optimizer
    nlopt::algorithm nlopt_algo;
    if (opt_name == "LN_COBYLA"){
        nlopt_algo = nlopt::LN_COBYLA;
    }
    else if (opt_name == "LN_BOBYQA"){
        nlopt_algo = nlopt::LN_BOBYQA;   // rounding errors
    }
    else if (opt_name == "LN_SBPLX"){
        nlopt_algo = nlopt::LN_SBPLX;    // failed to converge 2/10
    }
    else if (opt_name == "LN_NELDERMEAD"){
        nlopt_algo = nlopt::LN_NELDERMEAD;  
    }
    else if (opt_name == "LN_PRAXIS"){
        nlopt_algo = nlopt::LN_PRAXIS;  // just stops???
    }
    else {
        return 0;
    }

    nlopt::opt opt = nlopt::opt(nlopt_algo, params.all.size());
    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);
    opt.set_initial_step(steps);
    opt.set_ftol_abs(tolerance);

    std::vector<MOMAdata *> p_roots = get_roots(cells);

    opt.set_min_objective(target_func, &p_roots); // is type casted to void pointer (!)

    double ll_min;
    _save_ll = true;

    _file_log << "Optimization algorithm: " << opt.get_algorithm_name() << " Tolerance: " << tolerance << "\n";
    _iteration = 0;

    // actual minimization
    try{
        opt.optimize(parameter_state, ll_min);
    }

    catch(nlopt::roundoff_limited &e){
        _file_log << "(minimize_wrapper) WARNING: Log likelihood maximization is limited by rounding precision and was stopped. \
Although the tolerance criterium was not met, the last valid step is used for parameter estimation. (" << e.what() << ")\n";
    }

    catch(std::exception &e) {
        _file_log << "(minimize_wrapper) ERROR: Log likelihood optimization failed (" << e.what() << ")" << std::endl;
        throw;
    }
    _save_ll = false; // stop ll output

    _file_log << "Found maximum: log likelihood = " << std::setprecision(20) << -ll_min << std::setprecision(10) << "\n";
    params.set_final(parameter_state);
    _file_log << params << std::endl;
    return -ll_min;
}



// ================================================================================================ //
// ================================ in LOG parameter space ======================================== //
// ================================================================================================ //
double minimize_wrapper_log_params(double (*target_func)(const std::vector<double> &x, std::vector<double> &grad, void *p),
                        std::vector<MOMAdata> &cells,
                        Parameter_set &params, 
                        double tolerance, 
                        std::string opt_name){

    /* 
    * same as minimize_wrapper, but the parameter space is searched in log-space
    * the target function needs to be aware of the fact that the parameters are in log!
    */

    // set parameter space 
    std::vector<double> lower_bounds(params.all.size());
    std::vector<double> upper_bounds(params.all.size());
    std::vector<double> steps(params.all.size());

    std::vector<double> parameter_state(params.all.size());

    // everything needs to be in log
    // -------------------------------------------- //
    // -------------------------------------------- //
    for (size_t i=0; i<params.all.size(); ++i){
        if (params.all[i].minimized){
            parameter_state[i] = log(params.all[i].final);
        }
        else{
            parameter_state[i] = log(params.all[i].init);
        }
        if (params.all[i].fixed){
            steps[i] = 1.;
            lower_bounds[i] = log(params.all[i].init);
            upper_bounds[i] = log(params.all[i].init);
        } else{
            steps[i] = log(1.+params.all[i].step/params.all[i].init);
            lower_bounds[i] = log(params.all[i].lower);
            upper_bounds[i] = log(params.all[i].upper);
        }
    }
    // -------------------------------------------- /
    // -------------------------------------------- //

    // set up optimizer
    nlopt::algorithm nlopt_algo;
    if (opt_name == "LN_COBYLA"){
        nlopt_algo = nlopt::LN_COBYLA;
    }
    else if (opt_name == "LN_BOBYQA"){
        nlopt_algo = nlopt::LN_BOBYQA;   // rounding errors
    }
    else if (opt_name == "LN_SBPLX"){
        nlopt_algo = nlopt::LN_SBPLX;    // failed to converge 2/10
    }
    else if (opt_name == "LN_NELDERMEAD"){
        nlopt_algo = nlopt::LN_NELDERMEAD;  // failed to converge 2/10 (different ones)
    }
    else if (opt_name == "LN_PRAXIS"){
        nlopt_algo = nlopt::LN_PRAXIS;  // just stops???
    }
    else {
        return 0;
    }

    nlopt::opt opt = nlopt::opt(nlopt_algo, params.all.size());
    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);
    opt.set_initial_step(steps);
    opt.set_ftol_abs(tolerance);

    std::vector<MOMAdata *> p_roots = get_roots(cells);

    opt.set_min_objective(target_func, &p_roots); // is type casted to void pointer

    double ll_min;
    _save_ll = true;

    _file_log << "Optimization algorithm: " << opt.get_algorithm_name() << " Tolerance: " << tolerance << "\n";
    _iteration = 0;

   // actual minimization
    try{
        opt.optimize(parameter_state, ll_min);
    }

    catch(nlopt::roundoff_limited &e){
        _file_log << "(minimize_wrapper) WARNING: Log likelihood maximization is limited by rounding precision and was stopped. \
Although the tolerance criterium was not met, the last valid step is used for parameter estimation. (" << e.what() << ")\n";
    }

    catch(std::exception &e) {
        _file_log << "(minimize_wrapper) ERROR: Log likelihood optimization failed (" << e.what() << ")" << std::endl;
        throw;
    }
    _save_ll = false; // stop ll output

    _file_log << "Found maximum: log likelihood = " << std::setprecision(20) << -ll_min << std::setprecision(10) << "\n";

    // save final value for each parameter
    for(size_t i=0; i<parameter_state.size(); ++i){
        parameter_state[i] = exp(parameter_state[i]);
    }
    params.set_final(parameter_state);
    _file_log << params << std::endl;
    return -ll_min;
}

