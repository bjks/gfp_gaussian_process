#include "moma_input.h"
#include "mean_cov_model.h"
#include "Parameters.h"

#include <math.h>
#include <cmath>

/* 
* functions corresponding to backward part end with '_r'
*/

/* --------------------------------------------------------------------------
* FORWARD PREDICTION
* -------------------------------------------------------------------------- */

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

void posterior(Eigen::MatrixXd xgt, MOMAdata &cell, Eigen::Matrix2d S, Eigen::Matrix2d Si){
    // tested (i.e. same output as python functions)
    Eigen::MatrixXd K = cell.cov.block(0,0,2,4);
    cell.mean = cell.mean + K.transpose() * Si * xgt;
    cell.cov = cell.cov - K.transpose() * Si * K;
}

/* -------------------------------------------------------------------------- */
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
    * Recursive implementation that applies the function sc_prediction_forward to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    sc_prediction_forward(params_vec, *cell);
    prediction_forward_recr(params_vec, cell->daughter1);
    prediction_forward_recr(params_vec, cell->daughter2);
}

void prediction_forward(const std::vector<double> &params_vec, std::vector<MOMAdata> &cells){
    /* applies prediction to each cell going down the tree starting from all root cells */
    std::vector<MOMAdata *> p_roots = get_roots(cells);

    for(size_t i=0; i<p_roots.size(); ++i){
        prediction_forward_recr(params_vec,  p_roots[i]);
    }
}


/* --------------------------------------------------------------------------
* --------------------------------------------------------------------------
* BACKWARD PREDICTION
* --------------------------------------------------------------------------
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

void mean_cov_model_r(MOMAdata &cell, 
                    double t, double ml, 
                    double gl, double sl2, 
                    double mq, double gq, 
                    double sq2, double b){
    /* reverses the mean_cov_model function by switching the sign OU process params and beta */
    mean_cov_model(cell,t,-ml,-gl,sl2,-mq,-gq,sq2,-b);
}

void append_reversed_mean(MOMAdata &cell){
    /* append the "reverse" of the mean 
    mean ->     + + - - 
    in front of mean_backward variable in cell
    */
    Eigen::VectorXd temp_mean(4); 
    temp_mean << cell.mean;

    temp_mean(2) = - cell.mean(2);
    temp_mean(3) = - cell.mean(3);
    cell.mean_backward.insert(cell.mean_backward.begin(), temp_mean);
}

void append_reversed_cov(MOMAdata &cell){
    /* append the "reverse" of the vov 
    cov ->  + + - - 
            + + - - 
            - - + + 
            - - + + 
    in front of cov_backward variable in cell
    */
    Eigen::MatrixXd temp_cov(4,4);
    temp_cov << cell.cov;
    std::vector<std::vector<int>> entries   {{0,2},
                                            {0,3},
                                            {1,2},
                                            {1,3}};

    for(size_t k=0; k<entries.size(); ++k){
        temp_cov(entries[k][0], entries[k][1]) = - cell.cov(entries[k][0], entries[k][1]);
        temp_cov(entries[k][1], entries[k][0]) = - cell.cov(entries[k][1], entries[k][0]);
    }
    cell.cov_backward.insert(cell.cov_backward.begin(), temp_cov);
}


/* -------------------------------------------------------------------------- */

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

    for (long t=cell.time.size()-1; t>-1; --t ){
        xg(0) = cell.log_length(t) - cell.mean(0);
        xg(1) = cell.fp(t)         - cell.mean(1);

        S = cell.cov.block(0,0,2,2) + D;
        Si = S.inverse();

        posterior(xg, cell, S, Si); // updates mean/cov

        // save current mean/cov before (!) those are set for the next time point
        append_reversed_mean(cell);
        append_reversed_cov(cell);

        // previous time point:
        if (t>0) {
            mean_cov_model_r(cell, cell.time(t)-cell.time(t-1) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }
    }
}


void prediction_backward_recr(const std::vector<double> &params_vec, 
                    MOMAdata *cell){
    /*  
    * Recursive implementation that applies the function sc_prediction_backward to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    prediction_backward_recr(params_vec, cell->daughter1);
    prediction_backward_recr(params_vec, cell->daughter2);
    sc_prediction_backward(params_vec, *cell);
}

void prediction_backward(const std::vector<double> &params_vec, std::vector<MOMAdata> &cells){
    std::vector<MOMAdata *> p_roots = get_roots(cells);

    for(size_t i=0; i<p_roots.size(); ++i){
        prediction_backward_recr(params_vec,  p_roots[i]);
    }
}


void combine_predictions(std::vector<MOMAdata> &cells){
    /* combines foward and backward predictions by multiplying the gaussians of those predictions */
    Eigen::VectorXd temp_mean(4); 
    Eigen::MatrixXd temp_cov(4,4); 

    for(size_t i=0; i<cells.size();++i){
        for (size_t j=0; j<cells[i].time.size();++j ){
            temp_mean << cells[i].mean_forward[j];
            temp_cov << cells[i].cov_forward[j];

            multiply_gaussian(temp_mean, temp_cov, 
                                cells[i].mean_backward[j], cells[i].cov_backward[j]);
            cells[i].mean_prediction.push_back(temp_mean);
            cells[i].cov_prediction.push_back(temp_cov);
        }
    }
}

/* --------------------------------------------------------------------------
* OUTPUT
* -------------------------------------------------------------------------- */
std::string outfile_name_prediction(std::map<std::string, std::string> arguments, std::string suffix=""){
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]) + "_prediction" + suffix;
    return outfile + ".csv";
}


void output_upper_triangle(std::ofstream &file, Eigen::MatrixXd m){
    /* Comma seperated output of upper triangle of 4x4 Eigen::matrix*/
    std::vector<std::vector<int>> upper_triangle {{0,0}, {0,1}, {0,2}, {0,3},
                                                         {1,1}, {1,2}, {1,3},
                                                                {2,2}, {2,3},
                                                                       {3,3}};
                                                                       
    for (size_t k=0; k<upper_triangle.size(); ++k){
        if (k>0)
            file << ",";
        file << m(upper_triangle[k][0], upper_triangle[k][1]);
    }
}

void output_vector(std::ofstream &file, Eigen::VectorXd v){
    /* Comma seperated output of Eigen::vector */
    for (size_t k=0; k<v.size(); ++k){
        if (k>0)
            file << ",";
        file << v(k);
    }
}


void write_pretictions_to_file(const std::vector<MOMAdata> &cells, std::string outfile, 
                                Parameter_set& params, std::string direction="n"){        
    params.to_csv(outfile);

    std::ofstream file(outfile, std::ios_base::app);
    file << "\ncell_id,time,log_length,fp,";
    file    << "mean_x,mean_g,mean_l,mean_q,"
            <<   "cov_xx,cov_xg,cov_xl,cov_xq," 
                     << "cov_gg,cov_gl,cov_gq,"
                            << "cov_ll,cov_lq,"
                                   << "cov_qq\n";
    for(size_t i=0; i<cells.size();++i){
        for (size_t j=0; j<cells[i].mean_forward.size();++j ){
            file << cells[i].cell_id << "," << cells[i].time[j] << "," << cells[i].log_length[j] << "," << cells[i].fp[j] << ",";
            if(direction=="f"){
                output_vector(file, cells[i].mean_forward[j]);
                file << ",";  
                output_upper_triangle(file, cells[i].cov_forward[j]);
            } else if (direction=="b"){
                output_vector(file, cells[i].mean_backward[j]);
                file << ",";  
                output_upper_triangle(file, cells[i].cov_backward[j]);
            } else{
                output_vector(file, cells[i].mean_prediction[j]);
                file << ",";  
                output_upper_triangle(file, cells[i].cov_prediction[j]);
            }
            file << "\n"; 
        }
    }

    file.close();
}