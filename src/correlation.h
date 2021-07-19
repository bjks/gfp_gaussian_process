#include "predictions.h"
#include "Gaussians.h"

/*
* This header relies on the Gaussian classes in contrast to the rest of the code, 
* which just handles mean and covariances seperately
*/

std::vector<std::vector<Gaussian>> init_joint_matrix(std::vector<MOMAdata> &cells){
    /*
    * inits a vector of vectors of Gaussians, 
    * which size is chosen such that joints 
    * P(z_n+1, z_n) ... P(z_n+N, z_n) corresponding to the columns 
    * fit for the largest possible N 
    */
    size_t maxt = 0;
    for (size_t i=0; i<cells.size(); ++i){
        if (cells[i].time.size() > maxt){
            maxt = cells[i].time.size();
        }
    }
    std::vector<std::vector<Gaussian>> joint_matrix(maxt-1);
    return joint_matrix;
}

/* =========================================================== */
/* Calculation of a single joint over arbitrarily spaced point */
/* =========================================================== */

Gaussian include_measurement(Gaussian joint, 
                             Eigen::MatrixXd D, 
                             double x, double g){
    /* 
    * include the measurements x and g into the joint distrubtion, D is the diag matrix with the measurment variances
    * it follows the calculation of the posterior in the prediction part
    */

    Gaussian measurement_distr(joint.m.head(2), joint.C.block(0,0,2,2) + D);

    Eigen::VectorXd xg(2);

    xg(0) = x - measurement_distr.m(0);
    xg(1) = g - measurement_distr.m(1);

    Eigen::Matrix2d Si = measurement_distr.C.inverse();
    Eigen::MatrixXd K = joint.C.block(0,0,2,joint.C.cols());
    
    Gaussian new_gaussian(joint.m + K.transpose() * Si * xg, 
                            joint.C - K.transpose() * Si * K);

    // std::cout << "\nKT Si K: \n" <<  K.transpose() * Si * K;
    return new_gaussian;
}


Gaussian consecutive_joint(const std::vector<double> &params_vec, MOMAdata cell, size_t t){
    /* given a P(z_n | D_n) the joint P(z_n+1, z_n | D_n) is returned */
    cell.mean = cell.mean_forward[t];
    cell.cov = cell.cov_forward[t];

    // C_n (D_n)
    Eigen::VectorXd mean1 = cell.mean;
    Eigen::MatrixXd cov1 = cell.cov;

    // K_n+1,n (D_n), this is calculated from C_n and mu_n
    double dt = cell.time(t+1)-cell.time(t);
    Eigen::MatrixXd cross_cov = cross_cov_model(cell, dt, params_vec[0], 
                                                params_vec[1], params_vec[2], params_vec[3], 
                                                params_vec[4], params_vec[5], params_vec[6]); 

    // C_n+1 (D_n), calculated from C_n and mu_n
    mean_cov_model(cell, dt, params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]);
    Eigen::VectorXd mean2 = cell.mean;    
    Eigen::MatrixXd cov2 = cell.cov;

    // construct the joint over [z_n+1, z_n]
    Eigen::VectorXd mean_joint(8); 
    mean_joint << mean2, mean1;
    
    Eigen::MatrixXd cov_joint(8, 8); 
    cov_joint << vstack(hstack(cov2,                    cross_cov), 
                        hstack(cross_cov.transpose(),   cov1));

    Gaussian joint(mean_joint, cov_joint);
    return joint;
}


Affine_gaussian consecutive_conditional(const std::vector<double> &params_vec, MOMAdata cell, size_t t){
    /* given a P(z_n | D_n) the conditional P(z_n+1 | z_n D_n+1) is returned */
    cell.mean = cell.mean_forward[t];
    cell.cov = cell.cov_forward[t];

    /* ------------ prior joint ------------*/
    // C_n (D_n)
    Eigen::VectorXd mean1 = cell.mean;
    Eigen::MatrixXd cov1 = cell.cov;

    // K_n+1,n (D_n), this is calculated from C_n and mu_n
    double dt = cell.time(t+1)-cell.time(t);
    Eigen::MatrixXd cross_cov = cross_cov_model(cell, dt, params_vec[0], 
                                                params_vec[1], params_vec[2], params_vec[3], 
                                                params_vec[4], params_vec[5], params_vec[6]); 

    // C_n+1 (D_n), calculated from C_n and mu_n
    mean_cov_model(cell, dt, params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]);
    Eigen::VectorXd mean2 = cell.mean;    
    Eigen::MatrixXd cov2 = cell.cov;

    // construct the joint over [z_n, z_n+1]
    Eigen::VectorXd mean_joint(8); 
    mean_joint << mean1, mean2;
    
    Eigen::MatrixXd cov_joint(8, 8); 
    cov_joint << vstack(hstack(cov1,                    cross_cov.transpose()), 
                        hstack(cross_cov,   cov2));

    Gaussian joint(mean_joint, cov_joint);

    // get conditional P (z_n+1, z_n) writen as gaussian N(z_n+1| ... ) -> N(z_n| ... )
    Affine_gaussian conditional = seperate_gaussian(joint).conditional.transform();
    return conditional;
}


// Multiplication of conditional of P(z_n+2, y_n+2 | z_n+1, D_n+1) with marginal of P(z_n+1, z_n, D_n+1)
Eigen::VectorXd calc_x(Gaussian n1, Affine_gaussian n2){
    return n2.A * (n1.C + n2.A).inverse() * n1.m \
         + n1.C * (n1.C + n2.A).inverse() * n2.a;
}

Eigen::MatrixXd calc_X(Gaussian n1, Affine_gaussian n2){
    return n1.C * (n1.C + n2.A).inverse() * n2.F;
}

Eigen::MatrixXd calc_Y(Gaussian n1, Affine_gaussian n2){
    return n1.C * (n1.C + n2.A).inverse() * n2.A;
}


//Integration / "propagation"
Affine_gaussian propagation(Affine_gaussian a, Affine_gaussian x){
    Affine_gaussian n(  a.a + a.F * x.a, 
                        a.F*x.F, 
                        a.A + a.F * x.A * a.F.transpose());
    return n;
}


Gaussian next_joint(Gaussian joint, Affine_gaussian conditional){
    /* 
    * given the joint P(z_n+1, z_n | D_n+1) as well as the conditional P(z_n+2 | z_n+1, D_n+1) 
    * the posterior P(z_n+2, z_n | D_n+1) is returned
    */
    Seperated_gaussian joint_sep =seperate_gaussian(joint);
    // multiplication
    // first factor after multiplication
    Eigen::VectorXd x = calc_x(joint_sep.marginal, conditional);
    Eigen::MatrixXd X = calc_X(joint_sep.marginal, conditional);
    Eigen::MatrixXd Y = calc_Y(joint_sep.marginal, conditional);

    Affine_gaussian NX(x, X, Y);

    // second factor after multiplication
    Affine_gaussian G(  conditional.a, 
                        conditional.F, 
                        joint_sep.marginal.C + conditional.A);

    // transform and evaluate at mean mu_n+1 (D_n+1)
    Gaussian next_marginal = G.transform(joint_sep.marginal.m);

    // Integration over z_n+1
    Affine_gaussian next_conditional = propagation(joint_sep.conditional, NX);

    Seperated_gaussian new_joint(next_marginal, next_conditional);
    // Seperated_gaussian new_joint(marg, next_conditional);

    /* ------------ include next x and g ------------*/

    return new_joint.to_joint();
}


Gaussian incorporate_backward_prob(Seperated_gaussian joint, Eigen::VectorXd mean_backward, Eigen::MatrixXd cov_backward){
    Gaussian backward(mean_backward, cov_backward);
    Gaussian marginal = Gaussian::multiply(joint.marginal, backward);
    Seperated_gaussian new_joint(marginal, joint.conditional);
    return new_joint.to_joint();
}

/* =========================================================== */
/* "Main" of the joint_distribution calculation */
/* =========================================================== */

void calc_joint_distributions(const std::vector<double> &params_vec, MOMAdata &cell, 
                             size_t n1, size_t n2, std::vector<std::vector<Gaussian>> &joint_matrix){
    if (n1==n2){
        return;
    }
    Affine_gaussian conditional;
    Gaussian combined_joint; 

    Eigen::MatrixXd D(2,2);
    D <<  params_vec[7], 0, 0,  params_vec[8];

    /* P(z_n | D_n) -> P(z_n+1, z_n | D_n1) */
    Gaussian joint = consecutive_joint(params_vec, cell, n1); 
    /* include x and g -> P(z_n+1, z_n | D_n1+1) */
    joint = include_measurement(joint, D, cell.log_length(n1+1), cell.fp(n1+1)); 
    
    combined_joint = incorporate_backward_prob(seperate_gaussian(joint), cell.mean_backward[n1+1], cell.cov_backward[n1+1]);
    joint_matrix[0].push_back(combined_joint);

    for (size_t t=n1+1; t<n2-1; ++t){
        /* P(z_t+1 | z_t , D_t+1) */
        conditional = consecutive_conditional(params_vec, cell, t);

        /* move to P(z_t+1, z_n | D_t) */
        joint = next_joint(joint, conditional);
        /* move to P(z_t+1, z_n | D_t+1) */
        joint = include_measurement(joint, D, cell.log_length(t+1), cell.fp(t+1)); 
    
        // append the joint for (n, t+1)
        combined_joint = incorporate_backward_prob(seperate_gaussian(joint), cell.mean_backward[t+1], cell.cov_backward[t+1]);
        joint_matrix[t-n1].push_back(combined_joint);
    }
}


/* ======================================================== */
/* Looping over pairs of points for joint distr calculation */
/* ======================================================== */

void sc_joint_distributions(const std::vector<double> &params_vec, MOMAdata &cell,
                    std::vector<std::vector<Gaussian>> &joint_matrix){
    for (size_t t=0; t<cell.time.size()-1; ++t){
        calc_joint_distributions(params_vec, cell, t, cell.time.size(), joint_matrix);
    }
}


std::vector<std::vector<Gaussian>> collect_joint_distributions(const std::vector<double> &params_vec, std::vector<MOMAdata> &cells){
    /* 
    * has to be run after the prediction is run !!!
    */
    std::vector<std::vector<Gaussian>> joint_matrix = init_joint_matrix(cells);

    for (size_t i=0; i< cells.size(); ++i){
        sc_joint_distributions(params_vec, cells[i], joint_matrix);
    }

    return joint_matrix;
}


/* -------------------------------------------------------------------------- */
Eigen::MatrixXd covariance_from_joint(std::vector<Gaussian> joints, bool naive=false){
    int n = joints[0].C.cols();
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(n, n);
    for (size_t i=0; i<n; ++i){
        for (size_t j=0; j<n; ++j){

            // calculate covariance element <C_ij + m_i *m_j> - <m_i> * <m_j>
            double mimj = 0;
            double mi = 0;
            double mj = 0;

            double mimi = 0;
            double mjmj = 0;

            for(size_t k=0; k<joints.size(); ++k){
                if (naive){
                    mimj += joints[k].m(i)*joints[k].m(j);

                    mimi += pow(joints[k].m(i),2);
                    mjmj += pow(joints[k].m(j),2);
                }
                else{
                    mimj += joints[k].C(i,j) + joints[k].m(i)*joints[k].m(j); // sum over all C_ij + m_i * m_j
                    mimi += pow(joints[k].m(i),2) + joints[k].C(i,i);
                    mjmj += pow(joints[k].m(j),2) + joints[k].C(j,j);
                }
                mi += joints[k].m(i);
                mj += joints[k].m(j);
            }
            cov(i,j) = mimj/joints.size() -  mi/joints.size() * mj/joints.size();

            // Normalization
            cov(i,j) /= sqrt((mimi/joints.size() - pow(mi/joints.size(),2)) * \
                             (mjmj/joints.size() - pow(mj/joints.size(),2))); 
            
            // ...move to the next entry
        }
    }
    return cov;
}

std::vector<size_t> count_joints(std::vector<std::vector<Gaussian>> joint_matrix){
    std::vector<size_t> n;
    for(size_t dt=0; dt<joint_matrix.size(); ++dt){
        n.push_back(joint_matrix[dt].size());
    }
    return n;
}


std::vector<Eigen::MatrixXd> covariance_function(std::vector<std::vector<Gaussian>> joint_matrix){
    std::vector<Eigen::MatrixXd> covariances;
    for(size_t dt=0; dt<joint_matrix.size(); ++dt){
        covariances.push_back(covariance_from_joint(joint_matrix[dt]));
    }
    return covariances;
}


/* --------------------------------------------------------------------------
* OUTPUT
* -------------------------------------------------------------------------- */
std::string outfile_name_covariances(std::map<std::string, std::string> arguments, Parameter_set params){
    /* Filename for a prediction file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]) + outfile_param_code(params) + "_correlation";
    return outfile + ".csv";
}

void write_covariances_to_file(std::vector<Eigen::MatrixXd> covariances, double dt, std::vector<size_t> joint_number, std::string outfile, 
                                Parameter_set& params, const CSVconfig &config){    
    params.to_csv(outfile);
    Eigen::IOFormat CommaFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "\n");

    std::vector<std::string> string_entries = {"x(t+dt)", "g(t+dt)", "l(t+dt)", "q(t+dt)", "x(t)", "g(t)", "l(t)", "q(t)"};

    std::ofstream file(outfile, std::ios_base::app);
    file << "\ndt,joint_number";
    for(size_t m=0; m<string_entries.size(); ++m){
        for(size_t n=0; n<string_entries.size(); ++n){
            if (m<=n)
                file << ",R(" << string_entries[m] << string_entries[n] << ")";
        }
    }
    file << "\n";
    for(size_t i=0; i<covariances.size();++i){
        file << (i+1.)*dt << "," << joint_number[i] << ",";
        output_upper_triangle(file, covariances[i]);
        file << "\n";
    }
}