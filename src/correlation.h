#include "predictions.h"

// Matric utils
Eigen::MatrixXd hstack(Eigen::MatrixXd A, Eigen::MatrixXd B){
    /* horizontal stack of A and B */
    Eigen::MatrixXd C(A.rows(), A.cols()+B.cols());
    C << A, B;
    return C;
}

Eigen::MatrixXd vstack(Eigen::MatrixXd A, Eigen::MatrixXd B){
    /* vertical stack of A and B */
    Eigen::MatrixXd C(A.rows()+B.rows(), A.cols()); 
    C << A, B;
    return C;
}


/* ================================================= */
/* ================================================= */
/* ================================================= */

/* ============== Gaussian class ============== */
class Gaussian{
    /*  
    * A class describing a gaussian: N(x | m, C)
    */
public:
    Eigen::VectorXd m;
    Eigen::MatrixXd C;
    Gaussian() = default;
    Gaussian(Eigen::VectorXd m_, Eigen::MatrixXd C_) {
        m = m_;
        C = C_;
    }
    
    Gaussian multiply(Gaussian n1, Gaussian n2);
};

Gaussian Gaussian::multiply(Gaussian n1, Gaussian n2){
    /* multipication ignoring normalization (!) */
    m = n2.C*(n1.C + n2.C).inverse()*n1.m \
      + n1.C*(n1.C + n2.C).inverse()*n2.m;
    C = n1.C*(n1.C+n2.C).inverse()*n2.C;
    Gaussian n3(m, C);
    return n3;
}

/* ============== Affine_gaussian class ============== */
class Affine_gaussian{
    /*  
    * A class describing a gaussian with a transformed mean: N(y | a + F x, A)
    */
public:
    Eigen::VectorXd a;
    Eigen::MatrixXd F;
    Eigen::MatrixXd A;

    Affine_gaussian() = default;
    Affine_gaussian(Eigen::VectorXd a_, Eigen::MatrixXd F_, Eigen::MatrixXd A_){
        a = a_;
        F = F_;
        A = A_;
    }
    Affine_gaussian transform();
    Gaussian transform(Eigen::VectorXd y);
};

Affine_gaussian Affine_gaussian::transform(){
    /* 
    * transforms affine gaussian N(y|a+Fx,A)=N(a+Fx|y,A) -> N(x|a' + F' y, A')
    * looses normalization 
    */
    Eigen::VectorXd new_a = -F.inverse()*a;
    Eigen::MatrixXd new_F = F.inverse();
    Eigen::MatrixXd new_A = F.inverse() * A *F.inverse().transpose();
    Affine_gaussian n(new_a, new_F, new_A);
    return n;
}

Gaussian Affine_gaussian::transform(Eigen::VectorXd y){
    /* transforms affine gaussian N(y|a+Fx,A)=N(a+Fx|y,A) -> N(x|m,C) for a given y , looses normalization */
    Gaussian n(F.inverse()*(y-a), F.inverse()*A*F.inverse().transpose());
    return n;
}


/* ============== Seperated_gaussian class ============== */
class Seperated_gaussian{
    /*  
    * A class describing a joint gaussian as a product of 
    * mearginal and conditional in the form: N(x | m, C) N(y | a + F x, A)
    * It makes use of the Gaussian class and the 
    */
public:
    Gaussian marginal;
    Affine_gaussian conditional; 

    Seperated_gaussian(Gaussian marginal_, Affine_gaussian conditional_){
        marginal = marginal_;
        conditional = conditional_;
    }

    Gaussian to_joint(int n=8);
};

Gaussian Seperated_gaussian::to_joint(int n){
    /* rewrites the sperated gaussian (N(x | m, C) N(y | a + F x, A)) as N([x y]| m, C) whith matching m and C  */
    Eigen::VectorXd mean_joint(n); 
    mean_joint << marginal.m, conditional.a + conditional.F * marginal.m;

    Eigen::MatrixXd cov_joint(n, n); 
    cov_joint << vstack(hstack(marginal.C , 
                                marginal.C.transpose()*conditional.F.transpose()) , 
                        hstack(conditional.F * marginal.C, 
                                conditional.A + conditional.F * marginal.C.transpose() * conditional.F.transpose()));
    Gaussian joint(mean_joint, cov_joint);
    return joint;   
}


Seperated_gaussian seperate_gaussian(Gaussian joint, int n=4){
    Eigen::MatrixXd B = joint.C.bottomRightCorner(n,n);
    Eigen::MatrixXd K = joint.C.topRightCorner(n,n);
    Eigen::MatrixXd A = joint.C.topLeftCorner(n,n);

    Eigen::VectorXd a = joint.m.head(n);
    Eigen::VectorXd b = joint.m.tail(n);

    Gaussian n1(a, A);
    Affine_gaussian n2(b - K.transpose()*A.inverse()*a, 
                                K.transpose()*A.inverse(), 
                                B- K.transpose()*A.inverse()*K);
    Seperated_gaussian sep(n1, n2);
    return sep;
}


//Correlation
Eigen::MatrixXd correlation_from_covariance(Eigen::MatrixXd cov){
    Eigen::MatrixXd corr(cov.rows(), cov.cols());
    for (size_t i=0; i<cov.rows(); ++i){
        for (size_t j=0; j<cov.cols(); ++j){
            corr(i,j) = cov(i,j) / sqrt(cov(i,i) * cov(j,j) );
        }
    }
    return corr;
}
/* ================================================= */
/* ================================================= */
/* ================================================= */

Gaussian include_measurement(Gaussian distr, 
                             Gaussian measurement_distr, 
                             double x, double g){
    /* include the measurements x and g, this follows the calculation of the posterior in the prediction part*/
    Eigen::VectorXd xg(2);

    xg(0) = x - measurement_distr.m(0);
    xg(1) = g - measurement_distr.m(1);

    Eigen::Matrix2d Si = measurement_distr.C.inverse();

    Eigen::MatrixXd K = distr.C.block(0,0,2,distr.C.cols());
    // Eigen::MatrixXd K = hstack(distr.C.block(0,0,2,4), Eigen::MatrixXd::Zero(2, 4));
    // std::cout << K << "\n";


    Gaussian new_gaussian(distr.m + K.transpose() * Si * xg, 
                            distr.C - K.transpose() * Si * K);

    // std::cout << "\nKT Si K: \n" <<  K.transpose() * Si * K;
    return new_gaussian;
}


Gaussian consecutive_joint(const std::vector<double> &params_vec, MOMAdata cell, size_t t){
    /* given a P(z_n | D_n) the joint P(z_n+1, z_n | D_n+1) is returned */
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

    // construct the joint over [z_n+1, z_n]
    Eigen::VectorXd mean_joint(8); 
    mean_joint << mean2, mean1;
    
    Eigen::MatrixXd cov_joint(8, 8); 
    cov_joint << vstack(hstack(cov2,                    cross_cov), 
                        hstack(cross_cov.transpose(),   cov1));

    Gaussian joint(mean_joint, cov_joint);

    /* ------------ include next x and g ------------*/
    Eigen::MatrixXd D(2,2);
    D <<  params_vec[7], 0, 0,  params_vec[8];

    Gaussian measurement_distr(joint.m.head(2), joint.C.block(0,0,2,2) + D);
    // std::cout << "\ncov prior \n" << joint.C;

    joint = include_measurement(joint, measurement_distr, cell.log_length(t+1), cell.fp(t+1)); 
    // std::cout << "\ncov post:\n" << joint.C ;

    return joint;
}


Affine_gaussian consecutive_conditional(const std::vector<double> &params_vec, MOMAdata cell, size_t t){
    Gaussian joint = consecutive_joint(params_vec, cell, t);

    // write joint P(z_n+1, z_n | D_n) as a gaussian over the vector (z_n, z_n+1) instead of (z_n+1, z_n)
    Eigen::VectorXd mean_joint(8); 
    mean_joint << joint.m.tail(4), joint.m.head(4);

    // flip C_n+1 with C_n and transpose K_n+1,n and K_n,n+1
    Eigen::MatrixXd cov_joint(8, 8); 
    cov_joint << vstack(hstack(joint.C.block(4,4,4,4),              joint.C.block(0,4,4,4).transpose() ), 
                        hstack(joint.C.block(4,0,4,4).transpose(),  joint.C.block(0,0,4,4)));

    Gaussian flipped_joint(mean_joint, cov_joint);

    // get conditional P (z_n+1, z_n) writen as gaussian N(z_n+1| ... ) -> N(z_n| ... )
    Affine_gaussian conditional = seperate_gaussian(flipped_joint).conditional.transform();
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

    Seperated_gaussian next_joint(next_marginal, next_conditional);

    return next_joint.to_joint();
}


/*
* ===============================================================
* ===============================================================
* ===============================================================
*/
void calc_correlation_matrix(const std::vector<double> &params_vec, MOMAdata &cell, size_t n1, size_t n2){

    /* P(z_n+1, z_n | D_n+1) */
    Gaussian joint = consecutive_joint(params_vec, cell, n1); 
    cell.correlation.push_back(correlation_from_covariance(joint.C));
    // cell.correlation.push_back(joint.C);

   
    for (size_t t=n1+1; t<n2-1; ++t){
        /* P(z_n+2 | z_n+1 , D_n+2) */
        Affine_gaussian conditional = consecutive_conditional(params_vec, cell, t);

        // move to P(z_n+2, z_n | D_n+2)
        joint = next_joint(joint, conditional);
    
        // append the correlation matrix for (n, n+2)
        cell.correlation.push_back(correlation_from_covariance(joint.C));
        // cell.correlation.push_back(joint.C);

    }
}

/* -------------------------------------------------------------------------- */
void sc_correlation(const std::vector<double> &params_vec, MOMAdata &cell){
    
    calc_correlation_matrix(params_vec, cell, 0, cell.time.size());
}

void correlation_recr(const std::vector<double> &params_vec, 
                    MOMAdata *cell){
    /*  
    * 
    */
    if (cell == nullptr)
        return;
    sc_correlation(params_vec, *cell);
    correlation_recr(params_vec, cell->daughter1);
    correlation_recr(params_vec, cell->daughter2);
}

void correlation(const std::vector<double> &params_vec, std::vector<MOMAdata> &cells){
    /* has be run after the prediction is run */
    std::vector<MOMAdata *> p_roots = get_roots(cells);

    for(size_t i=0; i<p_roots.size(); ++i){
        correlation_recr(params_vec,  p_roots[i]);
    }
}

/* --------------------------------------------------------------------------
* OUTPUT
* -------------------------------------------------------------------------- */

std::string outfile_name_correlation(std::map<std::string, std::string> arguments, Parameter_set params){
    /* Filename for a prediction file */
    std::string outfile = out_dir(arguments);
    outfile += file_base(arguments["infile"]) + outfile_param_code(params) + "_correlation";
    return outfile + ".csv";
}

void write_correlations_to_file(const std::vector<MOMAdata> &cells, std::string outfile, 
                                Parameter_set& params, const CSVconfig &config){    
    params.to_csv(outfile);
    Eigen::IOFormat CommaFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "\n");

    std::vector<std::string> string_entries = {"x(t+dt)", "g(t+dt)", "l(t+dt)", "q(t+dt)", "x(t)", "g(t)", "l(t)", "q(t)"};

    std::ofstream file(outfile, std::ios_base::app);
    file << "\ncell_id,parent_id,time";
    for(size_t m=0; m<string_entries.size(); ++m){
        for(size_t n=0; n<string_entries.size(); ++n){
            if (m<=n)
                file << ",R(" << string_entries[m] << string_entries[n] << ")";
        }
    }
    file << "\n";
    for(size_t i=0; i<cells.size();++i){
        for (size_t j=0; j<cells[i].correlation.size();++j ){
            file << cells[i].cell_id << "," << cells[i].parent_id << "," 
                 << cells[i].time[j] * config.rescale_time << ","; 

            output_upper_triangle(file, cells[i].correlation[j]);
            file << "\n";
        }
    }
}