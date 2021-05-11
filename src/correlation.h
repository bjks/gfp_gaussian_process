#include "predictions.h"

// The terms that come from joint 
Eigen::VectorXd calc_a(Eigen::VectorXd mean_n, Eigen::VectorXd mean_n1, Eigen::MatrixXd K_nn1, Eigen::MatrixXd C_n1){
    return mean_n - K_nn1 * C_n1.inverse() * mean_n1; 
}

Eigen::MatrixXd calc_F(Eigen::MatrixXd K_nn1, Eigen::MatrixXd C_n1){
    return K_nn1 * C_n1.inverse();
}

Eigen::MatrixXd calc_A(Eigen::MatrixXd C_n, Eigen::MatrixXd K_nn1, Eigen::MatrixXd C_n1){
    return C_n - K_nn1 * C_n1.inverse() * K_nn1.transpose();
}

// The terms that come from including the next data point
Eigen::VectorXd calc_b(Eigen::MatrixXd C_n1, Eigen::MatrixXd M, Eigen::VectorXd mean_n1, Eigen::VectorXd y_n1){
    return C_n1*(M+C_n1).inverse()*y_n1 + M*(M+C_n1).inverse()*mean_n1;
}

Eigen::MatrixXd calc_B(Eigen::MatrixXd C_n1, Eigen::MatrixXd M){
    return M*(M+C_n1).inverse()*C_n1;
}

//
Eigen::MatrixXd calc_x(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::VectorXd a, Eigen::VectorXd b){
    return A*(B+A).inverse() * b + B*(B+A).inverse()*a;
}

// The terms that come from the integration step
Eigen::MatrixXd calc_X(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd F){
    return B*(B+A).inverse()*F;
}

Eigen::MatrixXd calc_Y(Eigen::MatrixXd A, Eigen::MatrixXd B){
    return B*(B+A).inverse()*A;
}


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

/* Gaussian classes */
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
};

class Seperated_gaussian{
    /*  
    * A class describing a joint gaussian as a product of mearginal and conditional in the form: N(x | m, C) N(y | a + F x, A)
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


Gaussian consecutive_joint(const std::vector<double> &params_vec, MOMAdata cell, double dt){
    /* given a P(z_n,D_n) in form of a cell with given mean and covariance, the joint P(z_n+1, z_n, D_n) is returned */

    // C_n (D_n)
    Eigen::VectorXd mean1 = cell.mean;
    Eigen::MatrixXd cov1 = cell.cov;

    // K_n+1,n (D_n), this is calculated from C_n and mu_n
    Eigen::MatrixXd cross_cov = cross_cov_model(cell, dt, params_vec[0], 
                                                params_vec[1], params_vec[2], params_vec[3], 
                                                params_vec[4], params_vec[5], params_vec[6]); 

    // C_n+1 (D_n), calculated from C_n and mu_n
    mean_cov_model(cell, dt, params_vec[0], 
                        params_vec[1], params_vec[2], params_vec[3], 
                        params_vec[4], params_vec[5], params_vec[6]);
    Eigen::VectorXd mean2 = cell.mean;    
    Eigen::MatrixXd cov2 = cell.cov;

    Eigen::VectorXd mean_joint(8); 
    mean_joint << mean2, mean1;
    
    Eigen::MatrixXd cov_joint(8, 8); 
    cov_joint << vstack(hstack(cov2, cross_cov), hstack(cross_cov.transpose(), cov1));

    Gaussian joint(mean_joint, cov_joint);
    return joint;
}


Gaussian include_measurement(Gaussian gaussian, double x, double g, Eigen::MatrixXd D){
    /* include the measurements x and g, this follows the calcul*/
    Eigen::VectorXd xg(2);

    xg(0) = x - gaussian.m(0);
    xg(1) = g - gaussian.m(1);

    Eigen::Matrix2d S = gaussian.C.block(0,0,2,2) + D;
    Eigen::Matrix2d Si = S.inverse();

    Eigen::MatrixXd K = gaussian.C.block(0,0,2,4);
    Gaussian new_gaussian(gaussian.m + K.transpose() * Si * xg, 
                            gaussian.C - K.transpose() * Si * K);
    return new_gaussian;
}


Gaussian consecutive_conditional(Gaussian joint, MOMAdata cell){
    Eigen::MatrixXd cov_inv = joint.C.inverse();
    Eigen::MatrixXd W = cov_inv.bottomRightCorner(4,4);
    cov_inv.bottomRightCorner(4,4) = W - cell.cov.inverse();

    Gaussian new_joint(joint.m, cov_inv.inverse());
    return new_joint;
}

Eigen::MatrixXd correlation_from_covariance(Eigen::MatrixXd cov){
    Eigen::MatrixXd corr(cov.rows(), cov.cols());
    for (size_t i=0; i<cov.rows(); ++i){
        for (size_t j=0; j<cov.cols(); ++j){
            corr(i,j) = cov(i,j) / sqrt(cov(i,i) * cov(j,j) );
        }
    }
    return corr;
}

void correlation(const std::vector<double> &params_vec, MOMAdata &cell, size_t t1, size_t t2){
    cell.mean = cell.mean_forward[t1];
    cell.cov = cell.cov_forward[t1];
    int t=t1;
    Eigen::MatrixXd D(2,2);
    D <<  params_vec[7], 0, 0,  params_vec[8];

    /* P(z_n+1, z_n, D_n) */
    Gaussian joint = consecutive_joint(params_vec, cell, cell.time(t+1)-cell.time(t)); 

    /* P(z_n+1, z_n, D_n) / P(z_n, D_n)  */
    Gaussian conditional = consecutive_conditional(joint, cell); 


    /* -> N(x | m, C) N(y | a + F x, A) */
    Seperated_gaussian joint_sep = seperate_gaussian(joint); // 
    /* include x and g */
    joint_sep.marginal = include_measurement(joint_sep.marginal, cell.log_length(t), cell.fp(t), D); 
    Gaussian posterior_joint = joint_sep.to_joint();

    Seperated_gaussian conditional_sep = seperate_gaussian(conditional); // 
    /* include x and g */
    conditional_sep.marginal = include_measurement(conditional_sep.marginal, cell.log_length(t), cell.fp(t), D); 

    // std::cout << joint.C << "\n";

    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
    std::cout << correlation_from_covariance(joint.C).format(CSVFormat) << "\n\n";

    // std::cout << conditional.C << "\n\n";

    // std::cout << conditional_sep.to_joint().C << "\n";

}

/* -------------------------------------------------------------------------- */
void sc_correlation(const std::vector<double> &params_vec, MOMAdata &cell){
    correlation(params_vec, cell, 0, 1);
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