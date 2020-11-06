#include "moma_input.h"
#include <math.h> 
#include <cmath>

#define _USE_MATH_DEFINES

/* -------------------------------------------------------------------------- */
void mean_cov_after_division(MOMAdata &cell, double var_dx, double var_dg){
    Eigen::MatrixXd F = Eigen::MatrixXd::Identity(4, 4);
    F(1,1) = 0.5;
    Eigen::Vector4d f(-log(2.), 0.0, 0.0, 0.0);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);
    D(0,0) = var_dx;
    D(1,1) = var_dg;

    cell.mean = F*cell.parent->mean + f;
    cell.cov = D + F*cell.parent->cov*F.transpose();
}


/* -------------------------------------------------------------------------- */
Eigen::MatrixXd rowwise_add(Eigen::MatrixXd m, Eigen::VectorXd v){
    Eigen::MatrixXd ones = Eigen::MatrixXd::Constant(1,m.cols(), 1);
    Eigen::MatrixXd m_new = m;

    for(int i=0; i < m.rows(); ++i ){
        m_new.row(i) += ones * v(i);
    } 
    return m_new;
}

double  log_likelihood(MOMAdata &cell, double var_dx, double var_dg){
    Eigen::Matrix2d D;
    D << var_dx, 0, 
              0, var_dg;
    Eigen::MatrixXd S = cell.cov.block(0,0,2,2) + D;
    Eigen::MatrixXd Si = S.inverse();

    Eigen::MatrixXd y(2, cell.fp.size());
    y << cell.length.transpose() , cell.fp.transpose();
    y = rowwise_add(y, cell.mean.head(2));
    Eigen::MatrixXd a = -0.5 * y.transpose() * Si * y;

    return log_likelihood = a(0) -0.5 * log(S.determinant()) - 2* log(2*M_PI);
}
/* -------------------------------------------------------------------------- */


void posterior(MOMAdata &cell, double var_dx, double var_dg){
    Eigen::Matrix2d D;
    D << var_dx, 0, 
              0, var_dg;

    Eigen::MatrixXd S = cell.cov.block(0,0,2,2) + D;
    Eigen::MatrixXd Si = S.inverse();

    Eigen::MatrixXd K = cell.cov.block(0,0,2,4);

    Eigen::MatrixXd y(2, cell.fp.size());
    y << cell.length.transpose() , cell.fp.transpose();
    y = rowwise_add(y, -cell.mean.head(2));

    cell.mean = cell.mean + K.transpose() * Si * y;
    cell.cov = cell.cov - K.transpose() * Si * K;
}

/* -------------------------------------------------------------------------- */
void sc_likelihood(const std::vector<double> &params_vec, 
                    MOMAdata &cell, 
                    double &total_likelihood){
    double likelihood = 0;
    if (cell.is_root()){
        /* use initial conditions for nm, nC */
    }
    else{
        /* update nm and nC that depend on mother cell */
        mean_cov_after_division(cell, 1, 1);
    }
    
    for (int i=0; i<params_vec.size();++i){
        total_likelihood += log_likelihood(cell, 1, 1);
        posterior(cell, 1, 1);

    // mean_cov_model
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

    // MOMAdata cell = *(MOMAdata *) c;
    double total_likelihood = 0;

    likelihood_recr(params_vec,  (MOMAdata *) c, total_likelihood);
    return total_likelihood;
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
