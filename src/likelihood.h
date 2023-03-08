#include "correlation_tree.h"

#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <iomanip>

int _iteration = 0;
int _print_level;
bool _save_ll;
std::ofstream _file_iteration;


Eigen::MatrixXd rowwise_add(Eigen::MatrixXd m, Eigen::VectorXd v){
    /*
    * Adds a constant to a row of a matrix, where the constants for each row is given as a vector 
    */
    Eigen::MatrixXd ones = Eigen::MatrixXd::Constant(1,m.cols(), 1);
    Eigen::MatrixXd m_new = m;

    for(int i=0; i < m.rows(); ++i ){
        m_new.row(i) += ones * v(i);
    } 
    return m_new;
}

double log_likelihood(Eigen::MatrixXd xgt, MOMAdata &cell, Eigen::MatrixXd S, Eigen::MatrixXd Si){
    /*
    * log likelihood
    */
    Eigen::MatrixXd a = -0.5 * xgt.transpose() * Si * xgt;
    return a(0) -0.5 * log(S.determinant()) - 2* log(2*M_PI);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void sc_likelihood(const std::vector<double> &params_vec, 
                    MOMAdata &cell, 
                    double &tl){
    /* Calculates the likelihood of a single cell (can be a root cell)
    * the params_vec contains paramters in the following (well defined) order:
    * {mean_lambda, gamma_lambda, var_lambda, mean_q, gamma_q, var_q, beta, var_x, var_g, var_dx, var_dg}
    * 0             1               2           3       4       5       6       7   8       9       10   
    */
    init_sc_distribution(cell, params_vec);

    Eigen::VectorXd xg(2);

    Eigen::MatrixXd D(2,2);

    Eigen::Matrix2d S;
    Eigen::Matrix2d Si;

    for (size_t t=0; t<cell.time.size(); ++t ){
        xg(0) = cell.log_length(t) - cell.mean(0);
        xg(1) = cell.fp(t)         - cell.mean(1);


        /* if chosen, noise in fp is scaled with sqrt of the fp content itself */
        if (cell.noise_model == "scaled"){
            D <<  params_vec[7], 0, 0, params_vec[8]*(cell.mean(1)+cell.fp_auto);
        }
        else {
            D <<  params_vec[7], 0, 0,  params_vec[8];
        }

        S = cell.cov.block(0,0,2,2) + D;
        Si = S.inverse();

        tl += log_likelihood(xg, cell, S, Si); // add to total_likelihood of entire tree  
           
        if (std::isnan(tl)){
            if (_save_ll){
                _file_iteration << _iteration + 1 << ",";
                for (size_t i=0; i<params_vec.size(); ++i){
                    _file_iteration << std::setprecision(20) << params_vec[i]  << ",";
                }
                _file_iteration << std::setprecision(30) << tl << std::setprecision(15) << "\n";
            }
            
            _file_log << "\n(sc_likelihood) ERROR: Log likelihood is Nan\n";
            _file_log << "____________________________________________\n";
            _file_log << "cell_id: " << cell.cell_id << ", at time " << cell.time(t) << "\n";
            _file_log << _iteration + 1 << ": ";
            for (size_t i=0; i<params_vec.size(); ++i){
                _file_log << params_vec[i]  << ", ";
            }

            _file_log << "ll=" <<  std::setprecision(10) << tl  << "\n";

            _file_iteration.close();

            throw std::domain_error("Likelihood is Nan");
        }

        posterior(xg, cell, S, Si); // updates mean/cov     

        if (t<cell.time.size()-1) {
            mean_cov_model(cell, cell.time(t+1)-cell.time(t) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }   
    }
}


/* --------------------------------------------------------------------------
* liklihood wrapping
* -------------------------------------------------------------------------- */

void likelihood_recr(const std::vector<double> &params_vec, 
                    MOMAdata *cell, 
                    double &tl){
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    sc_likelihood(params_vec, *cell, tl);
    likelihood_recr(params_vec, cell->daughter1, tl);
    likelihood_recr(params_vec, cell->daughter2, tl);
}


double total_likelihood(const std::vector<double> &params_vec, std::vector<double> &grad, void *c){
    /*
    * total_likelihood of cell trees, to be maximized
    */

    double tl = 0;
    // type cast the void vector back to vector of MOMAdata pointers
    std::vector<MOMAdata*> cells = *(std::vector<MOMAdata*> *) c;

    for(size_t i=0; i < cells.size(); ++i){
        if (cells[i]->is_root() ){ //this check should be redundant?!
            likelihood_recr(params_vec,  cells[i] , tl);
        }
    }
    ++ _iteration;

    /* Save state of iteration in outfile */
    if (_save_ll){
        _file_iteration << _iteration << ",";
        for (size_t i=0; i<params_vec.size(); ++i){
            _file_iteration << std::setprecision(20) << params_vec[i]  << ",";
        }
        _file_iteration << std::setprecision(30) << tl << std::setprecision(15) << "\n";
    }

    /* Print output dependend on set _print_level */
    if (_print_level>0){
            std::cout << _iteration << ": ";
            for (size_t i=0; i<params_vec.size(); ++i){
                std::cout << std::setprecision(20) << params_vec[i]  << ", ";
            }
            std::cout << "ll=" << std::setprecision(30) << tl << std::setprecision(15) << "\n";
    }
    return -tl;
}

double total_likelihood_log_params(const std::vector<double> &log_params_vec, std::vector<double> &grad, void *c){
    std::vector<double> params_vec(log_params_vec.size());
    for (size_t i=0; i<log_params_vec.size(); ++i){
        params_vec[i] = exp(log_params_vec[i]);
    }
    return total_likelihood(params_vec, grad, c);
}


double total_likelihood(const std::vector<double> &params_vec, std::vector<MOMAdata> &cells){
    std::vector<double> g;
    std::vector<MOMAdata *> p_roots = get_roots(cells);
    return -total_likelihood(params_vec, g, &p_roots);
}


/* --------------------------------------------------------------------------
* ERROR BARS
* -------------------------------------------------------------------------- */
Eigen::MatrixXd num_jacobian_ll(Parameter_set &params, 
                                std::vector<MOMAdata> &cells, 
                                double epsilon){
    /* numerical estimation of jacobian of log_likelihood */
    double lminus, lplus;
    std::vector<double> xminus, xplus;
    double h;
    int ii;
    std::vector<double> params_vec = params.get_final();

    std::vector<int> idx_non_fixed = params.non_fixed();
    Eigen::MatrixXd jacobian(1, idx_non_fixed.size());

    for(size_t i=0; i<idx_non_fixed.size(); ++i){
        ii = idx_non_fixed[i];

        h = std::max(params_vec[ii] * epsilon, 1e-13);

        xminus = params_vec;
        xminus[ii] = xminus[ii] - h;
        lminus = total_likelihood(xminus, cells);

        xplus = params_vec;
        xplus[ii] = xplus[ii] + h;
        lplus = total_likelihood(xplus, cells);
        jacobian(0,i) = (lplus - lminus)/(2.*h);
    }
    return jacobian;
}


Eigen::MatrixXd num_hessian_ll(double (*func)(const std::vector<double> &p, std::vector<MOMAdata> &c),
                                Parameter_set &params, std::vector<MOMAdata> &cells, double epsilon){
    /* Computes numerical hessian matrix of target function:
    Hij = [f(x + hi ei + hj ej) - f(x + hi ei - hj ej) - f(x - hi ei + hj ej) + f(x - hi ei - hj ej) ]/(4 hi hj) 
    */
    double lij, li_j, l_ij, l_i_j;
    std::vector<double> xij, xi_j, x_ij, x_i_j;
    double h1, h2;
    int ii, jj;

    std::vector<double> params_vec = params.get_final();

    std::vector<int> idx_non_fixed = params.non_fixed();
    Eigen::MatrixXd hessian(idx_non_fixed.size(), idx_non_fixed.size());

    // i,j are the indices of the matrix, less or equal in size than parameter number
    // ii, jj are the indices of the full paramter set
    for(size_t i=0; i<idx_non_fixed.size(); ++i){ 
        ii = idx_non_fixed[i]; // paramterer index
        for(size_t j=0; j<idx_non_fixed.size(); ++j){
            jj = idx_non_fixed[j];
            h1 = std::max(params_vec[ii] * epsilon, 1e-12);
            h2 = std::max(params_vec[jj] * epsilon, 1e-12);

            xij = params_vec;
            xij[ii] = xij[ii] + h1;
            xij[jj] = xij[jj] + h2;
            lij = func(xij, cells);
            
            xi_j = params_vec;
            xi_j[ii] = xi_j[ii] + h1;
            xi_j[jj] = xi_j[jj] - h2;
            li_j = func(xi_j, cells);

            x_ij = params_vec;
            x_ij[ii] = x_ij[ii] - h1;
            x_ij[jj] = x_ij[jj] + h2;
            l_ij = func(x_ij, cells);

            x_i_j = params_vec;
            x_i_j[ii] = x_i_j[ii] - h1;
            x_i_j[jj] = x_i_j[jj] - h2;
            l_i_j = func(x_i_j, cells);
            hessian(i,j) = (lij - li_j - l_ij + l_i_j)/ (4*h1*h2);
        }
    }
    return hessian;
}

std::vector<double> ll_error_bars(Parameter_set &params, std::vector<MOMAdata> &cells, double epsilon){
    /* returns the squared error bars on parameters */
    Eigen::MatrixXd hessian_inv = num_hessian_ll(total_likelihood, params, cells, epsilon).inverse();
    
    std::vector<double> error;
    for(int i=0; i<hessian_inv.rows(); ++i){
        error.push_back(-hessian_inv(i, i));
    }
    return error;
}

/* --------------------------------------------------------------------------
* OUTPUT
* -------------------------------------------------------------------------- */

void setup_outfile_likelihood(std::string outfile, Parameter_set params){
    params.to_csv(outfile);
    std::ofstream file(outfile,std::ios_base::app);
    file << "\nlog_likelihoods:\niteration,";
    for (size_t i=0; i<params.all.size(); ++i){
        file << params.all[i].name << ",";
    }
    file << "log_likelihood" <<"\n";
    file.close();
}

/* -------------------------------------------------------------------------- */

std::string outfile_name_minimization_process(std::map<std::string, std::string> arguments, 
                                                Parameter_set params, int segment){
    std::string outfile = out_dir(arguments);
    outfile += add_segment_to_filename(file_base(arguments["infile"]), segment) + outfile_param_code(params);
    return outfile + "_iterations.csv";
}

std::string outfile_name_minimization_final(std::map<std::string, std::string> arguments, 
                                            Parameter_set params, int segment){
    std::string outfile = out_dir(arguments);
    outfile +=  add_segment_to_filename(file_base(arguments["infile"]), segment) + outfile_param_code(params);
    return outfile + "_final.csv";
}

void save_final_likelihood(std::string outfile, 
                            std::vector<MOMAdata> const &cells, 
                            double ll_max, 
                            std::string min_algo, 
                            double tolerance, 
                            std::string search_space,
                            std::string noise_model,
                            std::string cell_division_model,
                            std::string version){
    std::ofstream file(outfile,std::ios_base::app);
    long ndata_points = count_data_points(cells);
    file << "\n";
    file << "n_data_points, " << ndata_points << "\n";
    file << "total_log_likelihoood," << std::setprecision(15) << ll_max << "\n";
    file << "norm_log_likelihoood," << std::setprecision(15) << ll_max/ndata_points << "\n";
    file << "optimization_algorithm," << min_algo << "\n";
    file << "tolerance," << tolerance << "\n";
    file << "search_space," << search_space << "\n";
    file << "noise_model," << noise_model << "\n";
    file << "cell_division_model," << cell_division_model << "\n";
    file << "version," << version << "\n";
    file.close();
}

void save_error_bars(std::string outfile, Parameter_set &params, std::vector<MOMAdata> &cells){
    /* calculates and saves the error bars on parameter estimates in outfile */
    std::ofstream file(outfile,std::ios_base::app);
    file << "\nerrors^2:";
    file << "\nepsilon";

    std::vector<double> eps {5e-2, 1e-2, 5e-3};
    std::vector<double> error;
    
    std::vector<int> idx_non_fixed = params.non_fixed();
    for(size_t i=0; i<idx_non_fixed.size(); ++i){ 
        file << "," << params.all[idx_non_fixed[i]].name;
    }
    file << "\n";

    for(size_t i=0; i<eps.size(); ++i ){
        error = ll_error_bars(params, cells, eps[i]);
        file << eps[i];
        for (size_t j=0; j<error.size(); ++j ){
            file << "," << error[j];
        }
        file << "\n";
    }
    file.close();
}

// ================================================================================================ //

std::string outfile_name_scan(std::map<std::string, std::string> arguments, std::string var, int segment){
    std::string outfile = out_dir(arguments);
    outfile += add_segment_to_filename(file_base(arguments["infile"]), segment) + "_scan_" + var;
    return outfile + ".csv";
}

// ================================================================================================ //
std::string outfile_name_parameter_file(std::map<std::string, std::string> arguments, 
                                        Parameter_set params, int segment){
    /* Filename for a parameter file */
    std::string outfile = out_dir(arguments);
    outfile +=  add_segment_to_filename(file_base(arguments["infile"]), segment) + outfile_param_code(params) + "_parameter_file";
    return outfile + ".txt";
}


void create_parameter_file(std::string outfile, Parameter_set& params){ 
    std::ofstream file(outfile);
    file << "# Generated parameter file with the final parameters that may be used for predictions\n";
    for (size_t i=0; i<params.all.size(); ++i){
        file  << params.all[i].name << " = " << params.all[i].final << "\n";    
    }
    file.close();
}
