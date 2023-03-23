#include "predictions.h"

/*
* This header relies on the Gaussian classes in contrast to the rest of the code, 
* which just handles mean and covariances seperately
*/


std::vector<std::tuple<std::string, double>> get_column_indices(std::vector <MOMAdata> &cells){
    /*
    */
    std::vector<std::tuple<std::string, double>> column_indices;
    for (size_t i=0; i<cells.size(); ++i){
        for (size_t t=0; t<cells[i].time.size(); ++t){
            // if (cells[i].is_root() && t==0){
            //     continue;
            // }
            std::tuple<std::string, double> indx(cells[i].cell_id, cells[i].time[t]);
            column_indices.push_back(indx);
        }
    }
    return column_indices;
}


void add_to_joint_vector(const std::vector<std::tuple<std::string, double>> &column_indices, 
                        std::vector<Gaussian> &joint_vector, Gaussian joint, 
                        std::string cell_id, double t){
    /*
    */
    for (size_t i=0; i<column_indices.size();++i){
        if (std::get<0>(column_indices[i])==cell_id && std::get<1>(column_indices[i])==t){
            joint_vector[i] = joint;
        }
    }
}


class Joint_vector{
    /*  
    * Vector of joints that is indexed by the tuple of (cell_id, time)
    */
public:
    std::vector<std::tuple<std::string, double>> column_indices;
    std::vector<Gaussian> joints;
    std::vector<bool> is_set;

    Joint_vector() = default;
    Joint_vector(std::vector<MOMAdata> &cells) {
        column_indices = get_column_indices(cells);
        Gaussian temp;
        for (size_t i=0; i<column_indices.size(); ++i){
            joints.push_back(temp);
            is_set.push_back(false);
        }
    }

    void add(Gaussian joint, std::string cell_id, double t);
    size_t count();
    void clear();
    void write(std::ostream &file, int n);
    void write_column_indices(std::ostream &file, int n);
};

void Joint_vector::add(Gaussian joint, std::string cell_id, double t){
    /*
    * Add a join to the joint vector based in the cell_id and the time t
    */
    for (size_t i=0; i<column_indices.size();++i){
        if (std::get<0>(column_indices[i])==cell_id && std::get<1>(column_indices[i])==t){
            joints[i] = joint;
            is_set[i] = true;
        }
    }
}

void Joint_vector::clear(){
    /*
    * Resets the joint vector, only keeping the column_indices 
    */
    Gaussian temp;
    for (size_t i=0; i<column_indices.size(); ++i){
        joints[i] = temp;
        is_set[i] = false;
    }
}

size_t Joint_vector::count(){
    size_t count = 0; 
    for (size_t i=0; i<column_indices.size(); ++i){
        count += is_set[i];
    }
    return count;
}

void Joint_vector::write(std::ostream &file, int n=44){
    /*
    * Write down all joints that are added to the joint vector. Empty columns for not set ones
    */
    for (size_t i=0; i<column_indices.size();++i){
        if (is_set[i]){
            file << ',';
            output_vector(file, joints[i].m);
            file << ',';
            output_upper_triangle(file, joints[i].C);
        }
        else{
            file << std::string(n, ',');
        }
    }
}

void Joint_vector::write_column_indices(std::ostream &file, int n=44){
    /*
    * Write the column incides as 2.1_33 as one line 
    */
    for (size_t i=0; i<column_indices.size();++i){
        file << std::get<0>(column_indices[i]) << '_' << std::get<1>(column_indices[i]);
        if (i==column_indices.size()-1){
            file << std::string(n-1, ',');
        }
        else{
            file << std::string(n, ',');
        }
    }
}

/* ==================================================================================================== */
/* Calculation of a single joint over arbitrarily spaced point */
/* ==================================================================================================== */

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

    return new_gaussian;
}

/* ==================================================================================================== */
/* Joint and conditional at cell division */
/* ==================================================================================================== */

Gaussian consecutive_joint_cell_division(const std::vector<double> &params_vec, MOMAdata cell, size_t t){
    /* given a P(z_n | D_n) the joint P(z_n+1, z_n | D_n) is returned where z_n is before and z_n+1 is after cell division */
    
    cell.mean = cell.mean_forward[t];
    cell.cov = cell.cov_forward[t];

    // C_n (D_n)
    Eigen::VectorXd mean1 = cell.mean;
    Eigen::MatrixXd cov1 = cell.cov;

    mean_cov_model(cell, cell.daughter1->time(0) - cell.time(t), 
                        params_vec[0], params_vec[1], params_vec[2], 
                        params_vec[3], params_vec[4], params_vec[5], 
                        params_vec[6]);

    double var_dx = params_vec[9];
    double var_dg = params_vec[10];

    Eigen::MatrixXd F = Eigen::MatrixXd::Identity(4, 4);
    F(1,1) = 0.5;
    Eigen::Vector4d f(-log(2.), 0.0, 0.0, 0.0);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);

    Gaussian joint;

    if (cell.cell_division_model == "binomial"){
        /* cov of t+1 */
        cell.cov(0,0) += var_dx;
        cell.cov(0,1) = cell.cov(1,0) = cell.mean(1)/2. * var_dx + cell.cov(0,1); 
        cell.cov(1,1) =   var_dx * (cell.mean(1)*cell.mean(1) + cell.cov(1,1))/2.\
                        + var_dg * cell.mean(1)/4. * (1-var_dx) \
                        + cell.cov(1,1)/4.;

        cell.cov(2,1) /= 2;
        cell.cov(1,2) /= 2;

        cell.cov(3,1) /= 2;
        cell.cov(1,3) /= 2;

        cell.mean = F*cell.mean + f;

        Eigen::MatrixXd cross_cov = cell.cov_forward[t];
        cross_cov(1, 0) /= 2.;
        cross_cov(1, 1) /= 2.;
        cross_cov(1, 2) /= 2.;
        cross_cov(1, 3) /= 2.;

        Eigen::VectorXd mean2 = cell.mean;    
        Eigen::MatrixXd cov2 = cell.cov;

        // construct the joint over [z_n+1, z_n]
        Eigen::VectorXd mean_joint(8); 
        mean_joint << mean2, mean1;
        
        Eigen::MatrixXd cov_joint(8, 8); 
        cov_joint << vstack(hstack(cov2,                    cross_cov), 
                            hstack(cross_cov.transpose(),   cov1));
        
        joint = Gaussian(mean_joint, cov_joint);
    }
    else{
        D(0,0) = var_dx;
        D(1,1) = var_dg;
        cell.mean = F*cell.mean + f;
        cell.cov = D + F * cell.cov * F.transpose();

        /* write joint as seperated gaussian of the conditional z_n+1| z_n 
        (the model itself is formulated like that) and the marginal over z_n */
        Affine_gaussian conditional(f, F, D);
        Gaussian marginal(cell.mean_forward[t],  cell.cov_forward[t]);

        /* calculate the joint ie the 8 dimensional gaussian N( [z_n, z_n+1]^T |..., ... )*/
        Seperated_gaussian joint_sep(marginal, conditional);
        Gaussian joint_yx = joint_sep.to_joint();
        /* ->  N( [z_n+1, z_n]^T |..., ... ) */
        joint = joint_yx.flip_xy();
    }
    return joint;
}


Affine_gaussian consecutive_conditional_cell_division(const std::vector<double> &params_vec, MOMAdata cell, size_t t){ 

    cell.mean = cell.mean_forward[t];
    cell.cov = cell.cov_forward[t];

    // C_n (D_n)
    Eigen::VectorXd mean1 = cell.mean;
    Eigen::MatrixXd cov1 = cell.cov;

    mean_cov_model(cell, cell.daughter1->time(0) - cell.time(t), 
                        params_vec[0], params_vec[1], params_vec[2], 
                        params_vec[3], params_vec[4], params_vec[5], 
                        params_vec[6]);

    double var_dx = params_vec[9];
    double var_dg = params_vec[10];

    Eigen::MatrixXd F = Eigen::MatrixXd::Identity(4, 4);
    F(1,1) = 0.5;
    Eigen::Vector4d f(-log(2.), 0.0, 0.0, 0.0);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);

    Gaussian joint;

    if (cell.cell_division_model == "binomial"){
        /* cov of t+1 */
        cell.cov(0,0) += var_dx;
        cell.cov(0,1) = cell.cov(1,0) = cell.mean(1)/2. * var_dx + cell.cov(0,1); 
        cell.cov(1,1) =   var_dx * (cell.mean(1)*cell.mean(1) + cell.cov(1,1))/2.\
                        + var_dg * cell.mean(1)/4. * (1-var_dx) \
                        + cell.cov(1,1)/4.;

        cell.cov(2,1) /= 2;
        cell.cov(1,2) /= 2;

        cell.cov(3,1) /= 2;
        cell.cov(1,3) /= 2;

        cell.mean = F*cell.mean + f;

        Eigen::MatrixXd cross_cov = cell.cov_forward[t];
        cross_cov(1, 0) /= 2.;
        cross_cov(1, 1) /= 2.;
        cross_cov(1, 2) /= 2.;
        cross_cov(1, 3) /= 2.;

        Eigen::VectorXd mean2 = cell.mean;    
        Eigen::MatrixXd cov2 = cell.cov;
        // construct the joint over [z_n, z_n+1]
        Eigen::VectorXd mean_joint(8); 
        mean_joint << mean1, mean2;
        
        Eigen::MatrixXd cov_joint(8, 8); 
        cov_joint << vstack(hstack(cov1,                    cross_cov.transpose()), 
                            hstack(cross_cov,   cov2));

        Gaussian joint(mean_joint, cov_joint);

        // get conditional P (z_n+1 | z_n) writen as gaussian N(z_n+1| ... ) -> N(z_n| ... )
        Affine_gaussian conditional = seperate_gaussian(joint).conditional.transform();
        return conditional;
    }
    else{
        /* create F, f, D */
        Eigen::MatrixXd F = Eigen::MatrixXd::Identity(4, 4);
        F(1,1) = 0.5;
        Eigen::Vector4d f(-log(2.), 0.0, 0.0, 0.0);
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(4, 4);

        D(0,0) = var_dx;
        D(1,1) = var_dg;

        /* write conditional as affine gaussian of z_n+1| z_n (the model itself is formulated like that)
        * writen as gaussian N(z_n+1| ... ) -> N(z_n| ... )
        */
        Affine_gaussian conditional(f, F, D);
        return conditional.transform();
    }
}

/* ==================================================================================================== */
/* Joint and conditional using the model when cells do not divide */
/* ==================================================================================================== */

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
    /* given a P(z_n | D_n) the conditional P(z_n+1 | z_n, D_n+1) is returned */
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

    // get conditional P (z_n+1 | z_n) writen as gaussian N(z_n+1| ... ) -> N(z_n| ... )
    Affine_gaussian conditional = seperate_gaussian(joint).conditional.transform();
    return conditional;
}

/* ==================================================================================================== */
/* Integration/iteration */
/* ==================================================================================================== */

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
    Seperated_gaussian joint_sep = seperate_gaussian(joint);
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

    return new_joint.to_joint();
}


Gaussian incorporate_backward_prob(Seperated_gaussian joint, 
                                    Eigen::VectorXd mean_backward, 
                                    Eigen::MatrixXd cov_backward, 
                                    std::vector<double> params_vec){
    /* Incorporates backward distribution to the joint distribtion and divides by the prior P(z_n) */

    /* define prior */
    Eigen::VectorXd mean_prior(4);
    mean_prior << 0, 0, params_vec[0], params_vec[3];

    Eigen::MatrixXd inv_cov_prior = Eigen::MatrixXd::Zero(4, 4);
    inv_cov_prior(2,2) = (2.*params_vec[1])/params_vec[2];
    inv_cov_prior(3,3) = (2.*params_vec[4])/params_vec[5];

    /* divide backward part by prior */
    Eigen::MatrixXd new_backward_cov = (cov_backward.inverse() - inv_cov_prior).inverse();
    Eigen::VectorXd new_backward_mean = new_backward_cov*(cov_backward.inverse()*mean_backward - \
                                                        inv_cov_prior*mean_prior);
    
    Gaussian backward(new_backward_mean, new_backward_cov);

    /* multiply forward and backward part */
    Gaussian marginal = Gaussian::multiply(joint.marginal, backward);
    Seperated_gaussian new_joint(marginal, joint.conditional);
    return new_joint.to_joint();
}

bool crosscovariance_is_small(Gaussian joint, double tolerance){
    for (size_t i=0; i<4; ++i){
        for (size_t j=4; j<8; ++j){
            if(abs( joint.C(i,j) / (joint.m(i)*joint.m(j)) ) > tolerance){
                return false;
            }
        }
    }
    return true;
}

/* =========================================================== */
/* "Main" of the joint_distribution calculation */
/* =========================================================== */

bool calc_joint_distributions(  const std::vector<std::vector<double>> &params_vecs, 
                                MOMAdata &cell, 
                                int n, 
                                Joint_vector &joint_vector,
                                double tolerance_joint){

    /* Starting from P(z_n+1, z_n | D_n) */

    Eigen::MatrixXd D(2,2);
    Gaussian combined_joint;
    Affine_gaussian conditional;
    size_t m;

    /* First iteration (m=1) includes y_n+1 and calculates P(z_n+2 | z_n+1, D_n+1 )*/
    for (m=1 ; n+m<cell.time.size(); ++m){

        std::vector<double> params_vec = params_vecs[cell.segment[n+m]];

        /* ------------ Posterior ------------ */ 
        /* include x and g -> P(z_n+1, z_n | D_n1+1) */        
        if (cell.noise_model == "scaled"){
            D <<  params_vec[7], 0, 0,  params_vec[8]*(cell.mean(1)+cell.fp_auto);
        }
        else {
            D <<  params_vec[7], 0, 0,  params_vec[8];
        }

        cell.joint = include_measurement(cell.joint, D, cell.log_length(n+m), cell.fp(n+m)); 

        combined_joint = incorporate_backward_prob(seperate_gaussian(cell.joint ), 
                                                    cell.mean_backward[n+m], 
                                                    cell.cov_backward[n+m], 
                                                    params_vec);
        
        if (crosscovariance_is_small(combined_joint, tolerance_joint)){
            return true;
        }
        joint_vector.add(combined_joint, cell.cell_id, cell.time[n+m]);
        
        if (n+m<cell.time.size()-1){
            /* ------------ New prior ------------ */
            conditional = consecutive_conditional(params_vec, cell, n+m); // P(z_n+2 | z_n+1, D_n+1 )
            cell.joint = next_joint(cell.joint, conditional);
        }

        else{
            /* -------- New prior for cell(s) after division ------- */
            if (cell.daughter1 != nullptr){
                conditional = consecutive_conditional_cell_division(params_vecs[cell.segment[cell.time.size()-1]], cell, n+m); // P(z_n+2 | z_n+1, D_n+1 )
                cell.joint = next_joint(cell.joint, conditional);
                
                cell.daughter1->joint = cell.joint;
            }
            if (cell.daughter2 != nullptr){
                cell.daughter2->joint = cell.joint;
            }
        }
    }
    return false;
}



/* ======================================================== */
/* Looping over pairs of points for joint distr calculation */
/* ======================================================== */

void joint_distributions_recr(  const std::vector<std::vector<double>> &params_vecs, 
                                MOMAdata *cell, 
                                int n, 
                                Joint_vector &joint_vector, 
                                double tolerance_joint, 
                                bool is_joint_at_division){
    /*  
    * Recursive implementation that applies the function calc_joint_distributions to every cell in the genealogy
    */
    if (cell == nullptr)
        return;
    if (!is_joint_at_division){
        bool stop = calc_joint_distributions(params_vecs, *cell, n, joint_vector, tolerance_joint);
        if (stop){
            return;
        }
    }
    joint_distributions_recr(params_vecs, cell->daughter1, -1, joint_vector, tolerance_joint, false);
    joint_distributions_recr(params_vecs, cell->daughter2, -1, joint_vector, tolerance_joint, false);
}


void sc_joint_distributions(const std::vector<std::vector<double>> &params_vecs, 
                            MOMAdata &cell,
                            std::ostream &file, 
                            Joint_vector &joint_vector,
                            double tolerance_joint){

    /* 
    * Calculates all joints starting from the time points in this cell. 
    * Loops over all time points in the cell and prints row of joint matrix to file.
    */

    for (size_t n=0; n<cell.time.size(); ++n){
        joint_vector.clear();

        // P(z_n+1, z_n | D_n) (using Theta_n)
        if (n<cell.time.size()-1){
            cell.joint = consecutive_joint(params_vecs[cell.segment[n]], cell, n);
            joint_distributions_recr(params_vecs, &cell, n, joint_vector, tolerance_joint, false);
        }
        else{
            if (cell.daughter1 != nullptr){
                cell.joint = consecutive_joint_cell_division(params_vecs[cell.segment[n]], cell, n);
                cell.daughter1->joint = cell.joint;
            }
            if (cell.daughter2 != nullptr){
                cell.daughter2->joint = cell.joint;
            }
            joint_distributions_recr(params_vecs, &cell, n, joint_vector, tolerance_joint, true);
        }

        /* ============= OUTPUT ============= */
        // if (cell.is_leaf() && n==cell.time.size()-1){
        //         continue;
        // }
        file << cell.cell_id << "," << cell.parent_id << "," << cell.time[n] ;
        joint_vector.write(file);
        file << "\n";
    }
}


void collect_joint_distributions(const std::vector<std::vector<double>> &params_vecs, 
                                std::vector<MOMAdata> &cells, 
                                std::ostream &out, 
                                double tolerance_joint){
    /* 
    * has to be run after the prediction is run!!!
    * setups the joint vector and loops over cells
    */

    Joint_vector joint_vector(cells);

    out << ",";
    joint_vector.write_column_indices(out);
    out << "\n";

    for (size_t i=0; i<cells.size(); ++i){
        _file_log << "-> cell: " << i << "\n";
        sc_joint_distributions(params_vecs, cells[i], out, joint_vector, tolerance_joint);
    }
}


// /* -------------------------------------------------------------------------- */
// Eigen::MatrixXd covariance_from_joint(std::vector<Gaussian> joints, size_t n, bool naive=false){
//     /* 
//     * Calculates the (normalized) covariance matrix R(zi, zj) from joints P(zi, zj | D) 
//     * If no joints are in vector, Nans are returned
//     */

//     Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(n, n);

//     // Handle empty vectors
//     if (joints.size()==0){
//         for (size_t i=0; i<n; ++i){
//             for (size_t j=0; j<n; ++j){
//                 cov(i,j) = std::numeric_limits<double>::quiet_NaN();
//             }
//         }
//         return cov;
//     }

//     // Compute C for non-empty vectors
//     for (size_t i=0; i<n; ++i){
//         for (size_t j=0; j<n; ++j){

//             // calculate covariance element <C_ij + m_i *m_j> - <m_i> * <m_j>
//             double mimj = 0;
//             double mi = 0;
//             double mj = 0;

//             double mimi = 0;
//             double mjmj = 0;

//             for(size_t k=0; k<joints.size(); ++k){
//                 if (naive){
//                     mimj += joints[k].m(i)*joints[k].m(j);

//                     mimi += pow(joints[k].m(i),2);
//                     mjmj += pow(joints[k].m(j),2);
//                 }
//                 else{
//                     mimj += joints[k].C(i,j) + joints[k].m(i)*joints[k].m(j); // sum over all C_ij + m_i * m_j
//                     mimi += pow(joints[k].m(i),2) + joints[k].C(i,i);
//                     mjmj += pow(joints[k].m(j),2) + joints[k].C(j,j);
//                 }
//                 mi += joints[k].m(i);
//                 mj += joints[k].m(j);
//             }
//             cov(i,j) = mimj/joints.size() -  mi/joints.size() * mj/joints.size();

//             // Normalization
//             cov(i,j) /= sqrt((mimi/joints.size() - pow(mi/joints.size(),2)) * \
//                              (mjmj/joints.size() - pow(mj/joints.size(),2))); 
            
//             // ...move to the next entry
//         }
//     }
//     return cov;
// }

// std::vector<Eigen::MatrixXd> covariance_function(std::vector<std::vector<Gaussian>> joint_matrix){
//     /* returns (normalized) correlation matrices from vector of vector joints */
//     std::vector<Eigen::MatrixXd> covariances;
//     for(size_t dt=0; dt<joint_matrix.size(); ++dt){
//         covariances.push_back(covariance_from_joint(joint_matrix[dt], 8));
//     }
//     return covariances;
// }

// std::vector<size_t> count_joints(std::vector<std::vector<Gaussian>> joint_matrix){
//     /* Counts the number of joints (number of pairs of data points) for each dt */
//     std::vector<size_t> n;
//     for(size_t dt=0; dt<joint_matrix.size(); ++dt){
//         n.push_back(joint_matrix[dt].size());
//     }
//     return n;
// }



/* --------------------------------------------------------------------------
* OUTPUT
* -------------------------------------------------------------------------- */
std::string outfile_name_covariances(std::map<std::string, std::string> arguments, std::vector<Parameter_set>& params_list){
    /* Filename for a covariance file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]);
    for (size_t i=0; i<params_list.size(); ++i){
        outfile += outfile_param_code(params_list[i]);
    }
    return outfile  + "_correlation" + ".csv";
}

/* Output correlation */

void write_covariances_to_file(std::vector<Eigen::MatrixXd> covariances, double dt, std::vector<size_t> joint_number, std::string outfile, 
                                std::vector<Parameter_set>& params_list, const CSVconfig &config){    
    // same formating as prediction files
    for(size_t i=0; i<params_list.size(); ++i){
    if (i==0)
        params_list[i].to_csv(outfile);
    else
        params_list[i].to_csv(outfile, std::ios_base::app);
}     

    Eigen::IOFormat CommaFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "\n");

    std::vector<std::string> string_entries = {"x(t+dt)", "g(t+dt)", "l(t+dt)", "q(t+dt)", "x(t)", "g(t)", "l(t)", "q(t)"};

    std::ofstream file(outfile, std::ios_base::app);
    file << "\ndt,joint_number";
    for(size_t m=0; m<string_entries.size(); ++m){
        for(size_t n=0; n<string_entries.size(); ++n){
            if (m<=n)
                file << ",R(" << string_entries[m] << string_entries[n] << ")"; // calling the correlation R to avoid confusion
        }
    }
    file << "\n";
    for(size_t i=0; i<covariances.size();++i){
        file << (i+1.)*dt << "," << joint_number[i] << ",";
        output_upper_triangle(file, covariances[i]);
        file << "\n";
    }
}

/* Output joints */

std::string outfile_name_joints(std::map<std::string, std::string> arguments, std::vector<Parameter_set>& params_list){
    /* Filename for a covariance file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]);
    for (size_t i=0; i<params_list.size(); ++i){
        outfile += outfile_param_code(params_list[i]);
    }
    return outfile  + "_joints" + ".csv";
}

void setup_outfile_joints(std::ostream &file, std::vector<Parameter_set>& params_list){
    for(size_t i=0; i<params_list.size(); ++i){
        params_list[i].to_csv(file);
    } 
    file << "\ncell_id,parent_id,time";
}