#include <iostream>
#include <vector>
#include <iomanip> 

#include <nlopt.hpp>

#include "Parameters.h"


double minimize_wrapper(double (*target_func)(const std::vector<double> &x, std::vector<double> &grad, void *p),
                        std::vector<MOMAdata> &cells,
                        Parameter_set &params, 
                        double tolerance, 
                        bool &found_min,
                        std::string opt_name){

    /* 
    * Wraps the minimization step using the nlopt library. 
    * The target function is the log_likelihood calculation returning -log_likelihood.
    * Returns the maximum (!) of the log_likelihood and sets the boolean found_min accoring 
    * to whether the optimization completed.
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
            steps[i] = 1.;
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
        nlopt_algo = nlopt::LN_NELDERMEAD;  // failed to converge 2/10 (different ones)
    }
    else if (opt_name == "LN_PRAXIS"){
        nlopt_algo = nlopt::LN_PRAXIS;  // just stops???
    }
    else {
        return 0;
    }

    nlopt::opt opt = nlopt::opt(nlopt_algo, params.all.size());         // robust, but rather slow near mnimimum

    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);
    opt.set_initial_step(steps);

    opt.set_ftol_abs(tolerance);

    std::vector<MOMAdata *> p_roots = get_roots(cells);

    opt.set_min_objective(target_func, &p_roots); // is type casted to void pointer

    double ll_min;
    _save_ll = true;

    std::cout << "Optimization algorithm: " << opt.get_algorithm_name() << " Tolerance: " << tolerance << "\n";
    // actual minimization
    try{
        opt.optimize(parameter_state, ll_min);
        std::cout << "Found maximum: log likelihood = " << std::setprecision(10) << -ll_min << "\n";

        if (std::isnan(ll_min)){
            std::cerr << "(minimize_wrapper) ERROR: Log likelihood optimization failed: Log likelihood is Nan:\
You may try more precise initial values or a smaller step size\n";
            return 0;
        }
        
        // save final value for each parameter
        params.set_final(parameter_state);
        std::cout << params << std::endl;
        found_min = true;
        _save_ll = false;
        return -ll_min;

    }
    catch(nlopt::roundoff_limited &e){
        std::cerr << "(minimize_wrapper) WARNING: Log likelihood maximization is limited by rounding precision and was stopped. \
Although the tolerance criterium was not met, the last valid step is used for parameter estimation. (" << e.what() << ")\n";
        // save final value for each parameter
        params.set_final(parameter_state);
        std::cout << params << std::endl;
        _save_ll = false;
        return -ll_min;
    }
    catch(std::exception &e) {
        std::cerr << "(minimize_wrapper) ERROR: Log likelihood optimization failed: " << e.what() << std::endl;
        _save_ll = false;
        return 0;
    }
}
