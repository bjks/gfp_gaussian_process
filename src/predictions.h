#include "likelihood.h"
/* 
* heavily re-uses functions from likelihood calculation,
* functions corresponding to backward part end with '_r'
*/

/* --------------------------------------------------------------------------
* FORWARD PREDICTION
* -------------------------------------------------------------------------- */

void sc_prediction_forward(const std::vector<double> &params_vec, 
                    MOMAdata &cell){
/* 
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
        posterior(xg, cell, S, Si); // updates mean/cov

        // save current mean/cov before (!) those are set for the next time point
        cell.mean_forward.push_back(cell.mean);
        cell.cov_forward.push_back(cell.cov);

        // next time point:
        if (t<cell.time.size()-1) {
            mean_cov_model(cell, cell.time(t+1)-cell.time(t) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }
    }
}


void prediction_forward_recr(const std::vector<double> &params_vec, 
                    MOMAdata *cell){
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    sc_prediction_forward(params_vec, *cell);
    prediction_forward_recr(params_vec, cell->daughter1);
    prediction_forward_recr(params_vec, cell->daughter2);
}

void prediction_forward(const std::vector<double> &params_vec, MOMAdata &cell){
    prediction_forward_recr(params_vec,  &cell);
}


/* --------------------------------------------------------------------------
* backward PREDICTION
* -------------------------------------------------------------------------- */


void multiply_gaussian(Eigen::VectorXd &m1, Eigen::MatrixXd &c1, Eigen::VectorXd m2, Eigen::MatrixXd c2){
    /* Multiply first gaussian with second one - inplace multiplication */

    Eigen::MatrixXd new_c1 = (c1.inverse() + c2.inverse()).inverse();

    m1 = new_c1 * c1.inverse() * m1  +  new_c1 * c2.inverse() * m2;
    c1 = new_c1;
}


void mean_cov_after_division_r(MOMAdata &cell, double var_dx, double var_dg){
    /*
    * mean and covariance matrix are updated as cell division occurs backward in time
    */
    Eigen::MatrixXd F = Eigen::MatrixXd::Identity(4, 4);
    F(1,1) = 2; // 2! 
    Eigen::Vector4d f(log(2.), 0.0, 0.0, 0.0); // plus sign
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);
    D(0,0) = var_dx;
    D(1,1) = var_dg;

    cell.mean = F*cell.daughter1->mean + f;
    cell.cov = D + F * cell.daughter1->cov * F.transpose();
    
    if (cell.daughter2 != nullptr){
        Eigen::Vector4d mean2 = F*cell.daughter2->mean + f;
        Eigen::MatrixXd cov2 = D + F * cell.daughter2->cov * F.transpose();

        multiply_gaussian(cell.mean, cell.cov, mean2, cov2);
    }
}

void reverse_mean_cov(MOMAdata &cell, 
                    double t, double ml, 
                    double gl, double sl2, 
                    double mq, double gq, 
                    double sq2, double b){
    mean_cov_model(cell, t, ml, 
                    gl,  sl2, 
                    mq,  gq, 
                    sq2,  b); 

    cell.mean(2) = - cell.mean(2);
    cell.mean(3) = - cell.mean(3);

    cell.cov(0,2) = - cell.cov(0,2); 
    cell.cov(0,3) = - cell.cov(0,3);
    cell.cov(1,2) = - cell.cov(1,2);
    cell.cov(1,3) = - cell.cov(1,3);
}


void mean_cov_model_r(MOMAdata &cell, 
                    double t, double ml, 
                    double gl, double sl2, 
                    double mq, double gq, 
                    double sq2, double b){
    mean_cov_model(cell,t,-ml,-gl,sl2,-mq,-gq,sq2,-b);
}


void sc_prediction_backward(const std::vector<double> &params_vec, 
                    MOMAdata &cell){
/* 
* the params_vec contains paramters in the following (well defined) order:
* {mean_lambda, gamma_lambda, var_lambda, mean_q, gamma_q, var_q, beta, var_x, var_g, var_dx, var_dg}
*/
    if (cell.is_leaf()){
        cell.mean = cell.mean_init;
        cell.cov = cell.cov_init;
    }
    else{
        // mean/cov is calculated from mother cell, does not depend on mean/cov of cell itself
        mean_cov_after_division_r(cell, params_vec[9], params_vec[10]);
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
        posterior(xg, cell, S, Si); // updates mean/cov

        // save current mean/cov before (!) those are set for the next time point
        cell.mean_forward.push_back(cell.mean);
        cell.cov_forward.push_back(cell.cov);

        // next time point:
        if (t<cell.time.size()-1) {
            mean_cov_model(cell, cell.time(t+1)-cell.time(t) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }
    }
}


void prediction_backward_recr(const std::vector<double> &params_vec, 
                    MOMAdata *cell){
    /*  
    * Recursive implementation that applies the function func to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    prediction_backward_recr(params_vec, cell->daughter1);
    prediction_backward_recr(params_vec, cell->daughter2);
    sc_prediction_backward(params_vec, *cell);

}

void prediction_backward(const std::vector<double> &params_vec, MOMAdata &cell){
    prediction_backward_recr(params_vec,  &cell);
}
