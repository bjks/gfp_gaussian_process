#include <iostream>
#include <vector>
#include <iomanip> 

#include <nlopt.hpp>

#include "Parameters.h"



double myvfunc(const std::vector<double> &parameters, std::vector<double> &grad, void *p){
    MOMAdata cell = *(MOMAdata *) p;
    std::cout << cell.cell_id << std::endl; 
    double sum=0;
    for (size_t i=0; i<parameters.size();++i){
        sum += (i+1)*parameters[i];
    } 
    return pow(sum, 2);
}

int minimize_wrapper(double (*target_func)(const std::vector<double> &x, std::vector<double> &grad, void *p),
                        std::vector<MOMAdata> &cells,
                        Parameter_set &params, 
                        double relative_tol){

    // set parameter space 
    std::vector<double> lower_bounds(params.all.size());
    std::vector<double> upper_bounds(params.all.size());
    std::vector<double> steps(params.all.size());

    std::vector<double> parameter_state(params.all.size());

    for (size_t i=0; i<params.all.size(); ++i){
        parameter_state[i] = params.all[i].init;
        if (params.all[i].fixed){
            steps[i] = 1; // will not be used anyway, but needs to be non-zero
            lower_bounds[i] = params.all[i].init;
            upper_bounds[i] = params.all[i].init;
        } else if (params.all[i].bound){
            steps[i] = params.all[i].step;
            lower_bounds[i] = params.all[i].lower;
            upper_bounds[i] = params.all[i].upper;
        } else {
            steps[i] = params.all[i].step;
            lower_bounds[i] = -HUGE_VAL;
            upper_bounds[i] = HUGE_VAL;
        }
    }

    // set up optimizer
    nlopt::opt opt(nlopt::LN_COBYLA, params.all.size());

    opt.set_lower_bounds(lower_bounds);
    opt.set_upper_bounds(upper_bounds);
    opt.set_initial_step(steps);
    opt.set_xtol_rel(relative_tol);

    std::vector<MOMAdata *> p_roots = get_roots(cells);

    opt.set_min_objective(target_func, &p_roots); // is type casted to void pointer

    double minf;
    // actual minimization
    try{
        opt.optimize(parameter_state, minf);
        std::cout << "Found minimum: log likelihood = " << std::setprecision(10) << minf << std::endl;

        // save final value for each parameter
        params.set_final(parameter_state);
        std::cout << params << std::endl;

    }
    catch(std::exception &e) {
        std::cout << "Nlopt failed: " << e.what() << std::endl;
        return -1;
    }
    return 1;
}
