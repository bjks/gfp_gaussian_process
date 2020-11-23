#include "moma_input.h"
#include "mean_cov_model.h"
#include "in_output.h"

#include <math.h>
#include <cmath>

#define _USE_MATH_DEFINES


/* -------------------------------------------------------------------------- */
void mean_cov_after_division(MOMAdata &cell, double var_dx, double var_dg){
    // tested (i.e. same output as python functions)
    /*
    * mean and covariance matrix are updated as cell division occurs, thus 
    * this function is applied to cells that do have parent cells
    */
    Eigen::MatrixXd F = Eigen::MatrixXd::Identity(4, 4);
    F(1,1) = 0.5;
    Eigen::Vector4d f(-log(2.), 0.0, 0.0, 0.0);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);
    D(0,0) = var_dx;
    D(1,1) = var_dg;

    cell.mean = F*cell.parent->mean + f;
    cell.cov = D + F * cell.parent->cov * F.transpose();
}


/* -------------------------------------------------------------------------- */
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

void posterior(Eigen::MatrixXd xgt, MOMAdata &cell, Eigen::Matrix2d S, Eigen::Matrix2d Si){
    // tested (i.e. same output as python functions)
    Eigen::MatrixXd K = cell.cov.block(0,0,2,4);
    cell.mean = cell.mean + K.transpose() * Si * xgt;
    cell.cov = cell.cov - K.transpose() * Si * K;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void sc_likelihood(const std::vector<double> &params_vec, 
                    MOMAdata &cell, 
                    double &total_likelihood){
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

    for (long t=0; t<cell.time.size(); ++t ){
        xg(0) = cell.log_length(t) - cell.mean(0);
        xg(1) = cell.fp(t)         - cell.mean(1);

        S = cell.cov.block(0,0,2,2) + D;
        Si = S.inverse();

        total_likelihood += log_likelihood(xg, cell, S, Si); // add to total_likelihood of entire tree     
        posterior(xg, cell, S, Si); // updates mean/cov        

        if (t<cell.time.size()-1) {
            mean_cov_model(cell, cell.time(t+1)-cell.time(t) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }
        if (std::isnan(total_likelihood)){
            break;
        }
    }
}


/* --------------------------------------------------------------------------
* liklihood wrapping
* -------------------------------------------------------------------------- */

void likelihood_recr(const std::vector<double> &params_vec, 
                    MOMAdata *cell, 
                    double &total_likelihood){
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    sc_likelihood(params_vec, *cell, total_likelihood);
    likelihood_recr(params_vec, cell->daughter1, total_likelihood);
    likelihood_recr(params_vec, cell->daughter2, total_likelihood);
}


double total_likelihood(const std::vector<double> &params_vec, std::vector<double> &grad, void *c){
    /*
    * total_likelihood of cell tree, to be maximized
    */

    double total_likelihood = 0;
    likelihood_recr(params_vec,  (MOMAdata *) c, total_likelihood);
    ++ _iteration;

    /* Save state of iteration in outfile */
    std::ofstream file(_outfile,std::ios_base::app);

    file << _iteration << ",";
    for (int i=0; i<params_vec.size(); ++i){
        file << params_vec[i]  << ",";
    }
    file << std::setprecision(10) << total_likelihood  << "\n";
    file.close();
    if (_print_level>0)
        std::cout << _iteration << ": " << total_likelihood << "\n";
    return -total_likelihood;
}

double total_likelihood(const std::vector<double> &params_vec, void *c){
    std::vector<double> g;
    return total_likelihood(params_vec, g,c);
}


/* ==========================================================================
* DEMO on how this work
* ========================================================================== */

// example function to illustrate how this works
void set_generation(const std::vector<double> &params_vec, MOMAdata &cell){
    if (cell.parent != nullptr){
        cell.generation = cell.parent->generation + 1;
    } else{
        cell.generation = 0;
    }
}


void print_generation_tree(const std::vector<double> &params_vec, MOMAdata &cell){

    // apply the set_generation function to root cell followed by the first generation etc...
    apply_down_tree(params_vec, cell, set_generation);                    

    /* output the genealogy with generation
        -> 20150630.5.4.124 generation: 0 -> 20150630.5.4.133 generation: 1 
        -> 20150630.5.4.140 generation: 2 -> 20150630.5.4.148 generation: 3
    */
    std::vector<std::vector<MOMAdata *> > cell_paths = get_genealogy_paths(cell);

    std::cout << std::endl;
    for (std::vector<MOMAdata *> path : cell_paths){
        for (MOMAdata * cell : path){
            std::cout << " -> " << cell->cell_id << " generation: " << cell->generation ;
        }
        std::cout << "\n\n";
    }

}
