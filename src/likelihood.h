#include "predictions.h"

#define _USE_MATH_DEFINES

int _iteration = 0;
int _print_level;
std::string _outfile_ll;


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
    // tested (i.e. same output as python functions)
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
*/
    if (cell.is_root()){
        cell.mean = cell.mean_init;
        cell.cov = cell.cov_init;
    }
    else{
        // mean/cov is calculated from mother cell, does not depend on mean/cov of cell itself
        mean_cov_after_division(cell, params_vec[9], params_vec[10]);
    }

    Eigen::VectorXd xg(2);

    Eigen::MatrixXd D(2,2);
    D <<  params_vec[7], 0, 0,  params_vec[8];

    Eigen::Matrix2d S;
    Eigen::Matrix2d Si;

    for (size_t t=0; t<cell.time.size(); ++t ){
        xg(0) = cell.log_length(t) - cell.mean(0);
        xg(1) = cell.fp(t)         - cell.mean(1);

        S = cell.cov.block(0,0,2,2) + D;
        Si = S.inverse();

        tl += log_likelihood(xg, cell, S, Si); // add to total_likelihood of entire tree     
        posterior(xg, cell, S, Si); // updates mean/cov        

        if (t<cell.time.size()-1) {
            mean_cov_model(cell, cell.time(t+1)-cell.time(t) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }
        if (std::isnan(tl)){
            break;
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
        if (cells[i]->is_root() ){
            likelihood_recr(params_vec,  cells[i] , tl);
        }
    }
    ++ _iteration;

    /* Save state of iteration in outfile */
    std::ofstream file(_outfile_ll, std::ios_base::app);

    file << _iteration << ",";
    for (size_t i=0; i<params_vec.size(); ++i){
        file << params_vec[i]  << ",";
    }
    file << std::setprecision(10) << tl  << "\n";
    file.close();

    /* Print output dependend on set _print_level */
    if (_print_level>0){
        if (_print_level==0)
            std::cout << _iteration << ": " << tl << "\n";
        if (_print_level>0){
            std::cout << _iteration << ": ";
                for (size_t i=0; i<params_vec.size(); ++i){
                    std::cout << params_vec[i]  << ", ";
                }
                std::cout << "ll=" << tl  << "\n";
        }
    }

    return -tl;
}

double total_likelihood(const std::vector<double> &params_vec, std::vector<MOMAdata> &cells){
    std::vector<double> g;
    std::vector<MOMAdata *> p_roots = get_roots(cells);
    return total_likelihood(params_vec, g, &p_roots);
}

/* --------------------------------------------------------------------------
* ERROR BARS
* -------------------------------------------------------------------------- */
Eigen::MatrixXd num_jacobian_ll(Parameter_set &params, std::vector<MOMAdata> &cells, double epsilon){
    double lminus, lplus;
    std::vector<double> xminus, xplus;
    double h;
    int ii;
    std::vector<double> params_vec = params.get_final();

    std::vector<int> idx_non_fixed = params.non_fixed();
    Eigen::MatrixXd jacobian(idx_non_fixed.size(), 1);

    for(size_t i=0; i<idx_non_fixed.size(); ++i){
        ii = idx_non_fixed[i];

        h = std::max(params_vec[ii] * epsilon, 1e-13);

        xminus = params_vec;
        xminus[ii] = xminus[ii] - h;
        lminus = total_likelihood(xminus, cells);

        xplus = params_vec;
        xplus[ii] = xplus[ii] + h;
        lplus = total_likelihood(xplus, cells);
        jacobian(i,0) = (lplus - lminus)/(2.*h);
    }
    return jacobian;
}

Eigen::MatrixXd num_jac_hessian_ll(Parameter_set &params, std::vector<MOMAdata> &cells, double epsilon){
    Eigen::MatrixXd jacobian = num_jacobian_ll(params, cells, epsilon);
    return jacobian* jacobian.transpose() ;
}


Eigen::MatrixXd num_hessian_ll(Parameter_set &params, std::vector<MOMAdata> &cells, double epsilon){
    /* Computes approx. of hessian matrix of log-likelihood 
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
            if (!params.all[i].fixed){
                h1 = std::max(params_vec[ii] * epsilon, 1e-13);
                h2 = std::max(params_vec[jj] * epsilon, 1e-13);

                xij = params_vec;
                xij[ii] = xij[ii] + h1;
                xij[jj] = xij[jj] + h2;
                lij = total_likelihood(xij, cells);
                
                xi_j = params_vec;
                xi_j[ii] = xi_j[ii] + h1;
                xi_j[jj] = xi_j[jj] - h2;
                li_j = total_likelihood(xi_j, cells);

                x_ij = params_vec;
                x_ij[ii] = x_ij[ii] - h1;
                x_ij[jj] = x_ij[jj] + h2;
                l_ij = total_likelihood(x_ij, cells);

                x_i_j = params_vec;
                x_i_j[ii] = x_i_j[ii] - h1;
                x_i_j[jj] = x_i_j[jj] - h2;
                l_i_j = total_likelihood(x_i_j, cells);
                hessian(i,j) = (lij - li_j - l_ij + l_i_j)/ (4*h1*h2);

            }
        }
    }
    std::cout << hessian << "\n";
    std::cout << hessian << "\n";

    return (0.5*hessian * 0.5*hessian.transpose());
    
}

std::vector<double> ll_error_bars(Parameter_set &params, std::vector<MOMAdata> &cells, double epsilon = 1e-8){
    Eigen::MatrixXd hessian_inv = num_jac_hessian_ll(params,cells, epsilon).inverse();
    std::cout << hessian_inv << "\n\n";

    std::vector<double> error;
    for(int i=0; i<hessian_inv.rows(); ++i){
        error.push_back(sqrt(hessian_inv(i, i)));
    }
    return error;
}


/* --------------------------------------------------------------------------
* OUTPUT
* -------------------------------------------------------------------------- */

void setup_outfile_likelihood(std::string outfile, Parameter_set params){
    params.to_csv(outfile);
    std::ofstream file(outfile,std::ios_base::app);
    file << "\nlikelihoods:\niteration,";
    for (size_t i=0; i<params.all.size(); ++i){
        file << params.all[i].name << ",";
    }
    file << "likelihood" <<"\n";
    file.close();
}

/* -------------------------------------------------------------------------- */
std::string outfile_name_minimization(std::map<std::string, std::string> arguments, Parameter_set params){
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]) + "_f";

    for(size_t i=0; i < params.all.size() ;++i){
        if (!params.all[i].bound && !params.all[i].fixed){
            outfile += std::to_string(i);
        }
    }
    outfile += "_b";
    for(size_t i=0; i < params.all.size(); ++i){
        if (params.all[i].bound){
            outfile += std::to_string(i);
        }
    }
    return outfile + ".csv";
}


std::string outfile_name_scan(std::map<std::string, std::string> arguments, std::string var){
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]) + "_scan_" + var;
    return outfile + ".csv";
}
