#include "moma_input.h"
#include "mean_cov_model.h"
#include "Parameters.h"
#include "utils.h"

#include <math.h>
#include <cmath>

std::string _noise_model;

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

void init_sc_distribution(MOMAdata &cell, 
                        const double &mean_lambda, 
                        const double &gamma_lambda, 
                        const double &var_lambda, 
                        const double &mean_q, 
                        const double &gamma_q, 
                        const double &var_q, 
                        const double &var_dx, 
                        const double &var_dg){
    if (cell.is_root()){
        // set means
        cell.mean(0) = cell.mean_init_forward(0);
        cell.mean(1) = cell.mean_init_forward(1);

        cell.mean(2) = mean_lambda;
        cell.mean(3) = mean_q;

        // set cov
        cell.cov(0,0) = cell.cov_init_forward(0,0);
        cell.cov(1,1) = cell.cov_init_forward(1,1);

        cell.cov(2,2) = var_lambda/(2.*gamma_lambda);
        cell.cov(3,3) = var_q/(2.*gamma_q);

    }
    else{
        mean_cov_after_division(cell, var_dx, var_dg);
    }
}

void posterior(Eigen::MatrixXd xgt, MOMAdata &cell, Eigen::Matrix2d S, Eigen::Matrix2d Si){
    // tested (i.e. same output as python functions)
    Eigen::MatrixXd K = cell.cov.block(0,0,2,4);
    cell.mean = cell.mean + K.transpose() * Si * xgt;
    cell.cov = cell.cov - K.transpose() * Si * K;
}


/* -------------------------------------------------------------------------- */
void sc_prediction_forward(const std::vector<std::vector<double>> &params_vecs, 
                    MOMAdata &cell){
/* 
* the params_vec contains paramters in the following (well defined) order:
* {mean_lambda, gamma_lambda, var_lambda, mean_q, gamma_q, var_q, beta, var_x, var_g, var_dx, var_dg}
*/

    int segment;

    if (cell.is_root()){
        segment = 0;
    }
    else{
        segment = cell.parent->segment[cell.parent->segment.size()-1];
    }

    init_sc_distribution(cell, 
                        params_vecs[segment][0], 
                        params_vecs[segment][1], 
                        params_vecs[segment][2], 
                        params_vecs[segment][3], 
                        params_vecs[segment][4], 
                        params_vecs[segment][5], 
                        params_vecs[segment][9], 
                        params_vecs[segment][10]);
    
    // if (cell.is_root()){
    //     cell.mean = cell.mean_init_forward;
    //     cell.cov = cell.cov_init_forward;
    // }
    //     int segment = cell.parent->segment[cell.parent->segment.size()-1];
    //     // mean/cov is calculated from mother cell, does not depend on mean/cov of cell itself
    //     mean_cov_after_division(cell, params_vecs[segment][9], params_vecs[segment][10]);
    // }

    Eigen::VectorXd xg(2);

    Eigen::MatrixXd D(2,2);

    Eigen::Matrix2d S;
    Eigen::Matrix2d Si;
    std::vector<double> params_vec;

    for (long t=0; t<cell.time.size(); ++t){
        params_vec = params_vecs[cell.segment[t]];

        xg(0) = cell.log_length(t) - cell.mean(0);
        xg(1) = cell.fp(t)         - cell.mean(1);

        /* if chosen, noise in fp is scaled with sqrt of the fp content itself */
        if (_noise_model == "scaled"){
            D <<  params_vec[7], 0, 0,  params_vec[8]*cell.mean(1);
        }
        else {
            D <<  params_vec[7], 0, 0,  params_vec[8];
        }

        S = cell.cov.block(0,0,2,2) + D;
        Si = S.inverse();
        posterior(xg, cell, S, Si); // updates mean/cov

        // save current mean/cov before (!) those they are changed again
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


void prediction_forward_recr(const std::vector<std::vector<double>> &params_vecs, 
                    MOMAdata *cell){
    /*  
    * Recursive implementation that applies the function sc_prediction_forward to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    sc_prediction_forward(params_vecs, *cell);
    prediction_forward_recr(params_vecs, cell->daughter1);
    prediction_forward_recr(params_vecs, cell->daughter2);
}

void prediction_forward(const std::vector<std::vector<double>> &params_vecs, std::vector<MOMAdata> &cells){
    /* applies prediction to each cell going down the tree starting from all root cells */
    std::vector<MOMAdata *> p_roots = get_roots(cells);

    for(size_t i=0; i<p_roots.size(); ++i){
        prediction_forward_recr(params_vecs,  p_roots[i]);
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


Eigen::VectorXd reverse_mean(Eigen::VectorXd mean){
    Eigen::VectorXd temp_mean(4); 

    temp_mean << mean;
    temp_mean(2) = - mean(2);
    temp_mean(3) = - mean(3);
    return temp_mean;
}

Eigen::MatrixXd reverse_cov(Eigen::MatrixXd cov){
    Eigen::MatrixXd temp_cov(4,4);

    temp_cov << cov;
    std::vector<std::vector<int>> entries   {{0,2},
                                            {0,3},
                                            {1,2},
                                            {1,3}};

    for(size_t k=0; k<entries.size(); ++k){
        temp_cov(entries[k][0], entries[k][1]) = - cov(entries[k][0], entries[k][1]);
        temp_cov(entries[k][1], entries[k][0]) = - cov(entries[k][1], entries[k][0]);
    }
    return temp_cov;
}


void posterior_r(Eigen::MatrixXd xgt, MOMAdata &cell, Eigen::Matrix2d S, Eigen::Matrix2d Si){
    cell.mean = reverse_mean(cell.mean);
    cell.cov = reverse_cov(cell.cov);

    Eigen::MatrixXd K = cell.cov.block(0,0,2,4);
    cell.mean = cell.mean + K.transpose() * Si * xgt;
    cell.cov = cell.cov - K.transpose() * Si * K;

    cell.mean = reverse_mean(cell.mean);
    cell.cov = reverse_cov(cell.cov);
}

void init_sc_distribution_r(MOMAdata &cell, 
                        const double &mean_lambda, 
                        const double &gamma_lambda, 
                        const double &var_lambda, 
                        const double &mean_q, 
                        const double &gamma_q, 
                        const double &var_q, 
                        const double &var_dx, 
                        const double &var_dg){
    if (cell.is_leaf()){
        // set means
        cell.mean(0) = cell.mean_init_backward(0);
        cell.mean(1) = cell.mean_init_backward(1);

        cell.mean(2) = - mean_lambda;
        cell.mean(3) = - mean_q;

        // set cov
        cell.cov(0,0) = cell.cov_init_backward(0,0);
        cell.cov(1,1) = cell.cov_init_backward(1,1);

        cell.cov(2,2) = var_lambda/(2.*gamma_lambda);
        cell.cov(3,3) = var_q/(2.*gamma_q);

    }
    else{
        mean_cov_after_division_r(cell, var_dx, var_dg);
    }
}


void mean_cov_model_r(MOMAdata &cell, 
                    double t, double ml, 
                    double gl, double sl2, 
                    double mq, double gq, 
                    double sq2, double b){
    /* reverses the mean_cov_model function by switching the sign OU process params and beta */
    mean_cov_model(cell, t, -ml, gl, sl2, -mq, gq, sq2, -b);
}


void append_reversed_mean(MOMAdata &cell){
    /* append the "reverse" of the mean 
    mean ->     + + - - 
    in front of mean_backward variable in cell
    */
    Eigen::VectorXd temp_mean(4); 
    temp_mean << reverse_mean(cell.mean);

    cell.mean_backward.insert(cell.mean_backward.begin(), temp_mean);
}

void append_reversed_cov(MOMAdata &cell){
    /* append the "reverse" of the cov 
    cov ->  + + - - 
            + + - - 
            - - + + 
            - - + + 
    in front of cov_backward variable in cell
    */
    Eigen::MatrixXd temp_cov(4,4);
    temp_cov << reverse_cov(cell.cov);

    cell.cov_backward.insert(cell.cov_backward.begin(), temp_cov);
}


/* -------------------------------------------------------------------------- */

void sc_prediction_backward(const std::vector<std::vector<double>> &params_vecs, 
                    MOMAdata &cell){
/* 
* the params_vec contains paramters in the following (well defined) order:
* {mean_lambda, gamma_lambda, var_lambda, mean_q, gamma_q, var_q, beta, var_x, var_g, var_dx, var_dg}
*/
    int segment;

    if (cell.is_leaf()){
        segment = 0;
    }
    else{
        segment = cell.segment[cell.segment.size()-1];
    }

    init_sc_distribution_r(cell, 
                        params_vecs[segment][0], 
                        params_vecs[segment][1], 
                        params_vecs[segment][2], 
                        params_vecs[segment][3], 
                        params_vecs[segment][4], 
                        params_vecs[segment][5], 
                        params_vecs[segment][9], 
                        params_vecs[segment][10]);

    // if (cell.is_leaf()){
    //     cell.mean = cell.mean_init_backward;
    //     cell.cov = cell.cov_init_backward;
    // }
    // else{
    //     int segment = cell.segment[cell.segment.size()-1];
    //     // mean/cov is calculated from mother cell, does not depend on mean/cov of cell itself
    //     mean_cov_after_division_r(cell, params_vecs[segment][9], params_vecs[segment][10]);
    // }

    Eigen::VectorXd xg(2);

    Eigen::MatrixXd D(2,2);

    Eigen::Matrix2d S;
    Eigen::Matrix2d Si;
    std::vector<double> params_vec;

    for (long t=cell.time.size()-1; t>-1; --t ){
        xg(0) = cell.log_length(t) - cell.mean(0);
        xg(1) = cell.fp(t)         - cell.mean(1);
        
        params_vec = params_vecs[cell.segment[t]];

        /* if chosen, noise in fp is scaled with sqrt of the fp content itself */
        if (_noise_model == "scaled"){
            D <<  params_vec[7], 0, 0,  params_vec[8]*cell.mean(1);
        }
        else {
            D <<  params_vec[7], 0, 0,  params_vec[8];
        }

        S = cell.cov.block(0,0,2,2) + D;
        Si = S.inverse();

        // save current mean/cov before (!) those are set for the next time point
        append_reversed_mean(cell);
        append_reversed_cov(cell);

        posterior(xg, cell, S, Si); // updates mean/cov

        // previous time point:
        if (t>0) {
            params_vec = params_vecs[cell.segment[t-1]];
            mean_cov_model_r(cell, cell.time(t)-cell.time(t-1) , params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]); // updates mean/cov
        }
    }
}


void prediction_backward_recr(const std::vector<std::vector<double>> &params_vecs, 
                    MOMAdata *cell){
    /*  
    * Recursive implementation that applies the function sc_prediction_backward to every cell in the genealogy
    * not meant to be called directly, see wrapper below
    */
    if (cell == nullptr)
        return;
    prediction_backward_recr(params_vecs, cell->daughter1);
    prediction_backward_recr(params_vecs, cell->daughter2);
    sc_prediction_backward(params_vecs, *cell);
}

void prediction_backward(const std::vector<std::vector<double>> &params_vecs, std::vector<MOMAdata> &cells){
    std::vector<MOMAdata *> p_roots = get_roots(cells);

    for(size_t i=0; i<p_roots.size(); ++i){
        prediction_backward_recr(params_vecs,  p_roots[i]);
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

std::string outfile_param_code(const Parameter_set params){
    std::string code = "_f";
     for(size_t i=0; i < params.all.size() ;++i){
        if (!params.all[i].bound && !params.all[i].fixed){
            code += std::to_string(i);
        }
    }
    code += "_b";
    for(size_t i=0; i < params.all.size(); ++i){
        if (params.all[i].bound){
            code += std::to_string(i);
        }
    }
    return code;
}


std::string outfile_name_prediction(std::map<std::string, std::string> arguments, 
                                    std::vector<Parameter_set> params_list, std::string suffix=""){
    /* Filename for a prediction file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]);
    for (size_t i=0; i<params_list.size(); ++i){
        outfile += outfile_param_code(params_list[i]);
    }
    return outfile  + "_prediction" + suffix + ".csv";
}

std::string outfile_name_prediction_segments(std::map<std::string, std::string> arguments, std::string suffix=""){
    /* Filename for a prediction file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]) + + "segments_prediction" + suffix;
    return outfile + ".csv";
}

void output_upper_triangle(std::ostream &file, const Eigen::MatrixXd mat){
    /* Comma seperated output of upper triangle of Eigen::matrix*/
    for(size_t m=0; m<mat.rows();++m){
        for(size_t n=0; n<mat.cols();++n){
            if (m<=n){
                if (n+m>0)
                    file << ",";
                file << mat(m, n);
            }
        }    
    }                     
}

void output_vector(std::ostream &file, const Eigen::VectorXd v){
    /* Comma seperated output of Eigen::vector */
    for (size_t k=0; k<v.size(); ++k){
        if (k>0)
            file << ",";
        file << v(k);
    }
}


void write_predictions_to_file(const std::vector<MOMAdata> &cells, std::string outfile, 
                                std::vector<Parameter_set> &params_list, std::string direction="n"){   
    for(size_t i=0; i<params_list.size(); ++i){
        if (i==0)
            params_list[i].to_csv(outfile);
        else
            params_list[i].to_csv(outfile, std::ios_base::app);
    }     

    std::ofstream file(outfile, std::ios_base::app);
    file << "\ncell_id,parent_id,time,log_length,fp,";
    file    << "mean_x,mean_g,mean_l,mean_q,"
            <<   "cov_xx,cov_xg,cov_xl,cov_xq," 
                     << "cov_gg,cov_gl,cov_gq,"
                            << "cov_ll,cov_lq,"
                                   << "cov_qq\n";
    for(size_t i=0; i<cells.size();++i){
        for (size_t j=0; j<cells[i].mean_forward.size();++j ){
            file << cells[i].cell_id << "," << cells[i].parent_id << "," 
                 << cells[i].time[j] << "," 
                 << cells[i].log_length[j] << "," << cells[i].fp[j] << ",";
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