#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Core>

// Matrix utils
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
    
    static Gaussian multiply(Gaussian n1, Gaussian n2);
    Gaussian flip_xy(int n);
};

Gaussian Gaussian::multiply(Gaussian n1, Gaussian n2){
    /* multipication ignoring normalization (!) */
    Eigen::VectorXd m = n2.C*(n1.C + n2.C).inverse()*n1.m \
      + n1.C*(n1.C + n2.C).inverse()*n2.m;
    Eigen::MatrixXd C = n1.C*(n1.C+n2.C).inverse()*n2.C;
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
    Eigen::MatrixXd new_A = F.inverse() * A * F.inverse().transpose();
    Affine_gaussian n(new_a, new_F, new_A);
    return n;
}

Gaussian Affine_gaussian::transform(Eigen::VectorXd y){
    /* transforms affine gaussian N(y|a+Fx,A)=N(a+Fx|y,A) -> N(x|m,C) for a given y , looses normalization */
    Gaussian n(F.inverse()*(y-a), F.inverse() * A *F.inverse().transpose());
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
    /* 
    * rewrites the sperated gaussians (N(x | m, C) N(y | a + F x, A)) as N([x y]| m, C) whith matching m and C  
    * "Inverse" of seperate_gaussian()
    */
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
    /* 
    * rewrites the joint  N([x y]| m, C) as sperated gaussians (N(x | m, C) N(y | a + F x, A)) whith matching m, C, a, F, A  
    * "Inverse" of to_joint() (and thus not part of Seperated_gaussian classe)
    */
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

Gaussian Gaussian::flip_xy(int n=4){
    /* flips x and y of 2n dimensional gaussian ie N([x,y]) -> N([y,x]), where x and y are n dimensional */
    Eigen::VectorXd mean(n*2); 
    mean << m.tail(n),  m.head(n);
    
    Eigen::MatrixXd cov(n*2, n*2); 
    cov << vstack(hstack(C.block(n,n,n,n),                C.block(0,n,n,n).transpose()), 
                  hstack(C.block(n,0,n,n).transpose(),    C.block(0,0,n,n)));

    Gaussian n_new(mean, cov);
    return n_new;
}