#include <cmath>
#include "Faddeeva.hh"

#define _USE_MATH_DEFINES

// ======================================================================================================== //       
// The following intergals are calculated via Dawson functions, to avoid numerical overflows

double zerotauint(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 exp[a*s**2+b*s+c]ds rewritten with Dawson functions 
    double x = 2. *(\
            - exp(a*pow(t0,2) + b*t0 +c) *\
             Faddeeva::Dawson((b + 2.*a*t0)/(2.*sqrt(a)) ) \
            + exp(a*pow(t1,2) + b*t1 +c) *\
            (Faddeeva::Dawson((b + 2.*a*t1)/(2.*sqrt(a)))) );

    // if (std::isnan(x)){
    //     std::cout   << "INF-WARNING: zerotauint (a,b,c,t0,t1) " << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "  ";
    // }
    return x /(2.*sqrt(a));
}

double onetauint(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 s*exp[a*s**2+b*s+c]ds rewritten with Dawson functions 
    double x= (-2.*sqrt(a) * exp(c) * (exp(t0*(b + a*t0)) - exp(t1*(b + a*t1))) +\
           b*2.* \
           (exp(-pow(b,2.)/(4.*a) + c + pow((b + 2.*a*t0)/(2.*sqrt(a)), 2.)) * Faddeeva::Dawson((b + 2.*a*t0)/(2.*sqrt(a))) \
           -exp(-pow(b,2.)/(4.*a) + c + pow((b + 2.*a*t1)/(2.*sqrt(a)), 2.)) * Faddeeva::Dawson((b + 2.*a*t1)/(2.*sqrt(a)))));

    // if (std::isnan(x)){
    //     std::cout   << "INF-WARNING: onetauint (a,b,c,t0,t1) " << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "\n";
    // }    
    return x/(4.*pow(a,1.5));
}

double twotauint(double a, double b, double c, double t_1, double t_0=0){
    //int_t_0^t_1 s**2*exp[a*s**2+b*s+c]ds rewritten with Dawson functions 
    double x = (2.*sqrt(a)*exp(c)*(exp(t_0*(a*t_0 + b))*(b - 2.*a*t_0) - exp(t_1*(a*t_1 + b))*(b - 2.*a*t_1))  +\
            (exp(-pow(b,2.)/(4.*a) + c + pow((b + 2.*a*t_0)/(2.*sqrt(a)),2.)) * \
            (2.*a - pow(b,2.))*2.*Faddeeva::Dawson((b + 2.*a*t_0)/(2.*sqrt(a))) +\
            exp(-pow(b,2)/(4.*a) + c + pow((b + 2.*a*t_1)/(2.*sqrt(a)),2.)) * \
            (-2.*a + pow(b,2.))*2.*Faddeeva::Dawson((b + 2.*a*t_1)/(2.*sqrt(a)))));

    // if (std::isnan(x)){
    //     std::cout  << "INF-WARNING: twotauint (a,b,c,t_0,t_1) " << a << " "<< b << " "<< c << " "<<  t_0 <<  " "<<  t_1 << "\n";
    //     std::cout  << exp(-pow(b,2)/(4.*a) + c + pow((b + 2*a*t_1)/(2.*sqrt(a)),2)) << "\n";
    // }
    return x/(8.*pow(a,2.5));
}

double treetauint(double a, double b, double c, double t_1, double t_0=0){
    //int_t_0^t_1 s**3*exp[a*s**2+b*s+c]ds rewritten with Dawson functions 
    double x = (-2.*sqrt(a)*exp(c)*\
            (pow(b,2.)*(exp(t_0*(b + a*t_0)) - exp(t_1*(b + a*t_1))) - \
            2.*a*exp(t_0*(b + a*t_0))*(2.+b*t_0) + 2*a*exp(t_1*(b + a*t_1))*(2 + b*t_1) +\
            4*pow(a,2.)*(exp(t_0*(b + a*t_0))*pow(t_0,2.) - exp(t_1*(b + a*t_1))*pow(t_1,2))))  + \
            exp(-pow(b,2)/(4.*a) + c + pow((b + 2.*a*t_0)/(2.*sqrt(a)), 2.))*b*(-6.*a + pow(b,2.))*2.* \
            Faddeeva::Dawson((b + 2.*a*t_0)/(2.*sqrt(a))) -\
            exp(-pow(b,2.)/(4.*a) + c + pow((b + 2.*a*t_1)/(2.*sqrt(a)), 2))*b*(-6*a + pow(b,2.))*2.* \
            Faddeeva::Dawson((b + 2.*a*t_1)/(2.*sqrt(a)));

    // if (std::isnan(x)){
    //     std::cout << "INF-WARNING: treetauint (a,b,c,t_0,t_1) " << a << " "<< b << " "<< c << " "<<  t_0 << " "<<  t_1 << "\n";
    //     std::cout << exp(-pow(b,2)/(4.*a) + c + pow((b + 2*a*t_1)/(2.*sqrt(a)), 2)) << "\n";
    // }
    return x/(16.*pow(a,3.5));
}

// ======================================================================================================== //
// MEAN ELEMENTS
// ======================================================================================================== //       

double mean_x(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return bx+ml*t+(bl-ml)*(1-exp(-gl*t))/gl;
}

double mean_g(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return bg/exp(b*t)+Clq*onetauint(Cll/2.,b+bl+Cxl-gq,bx+Cxx/2.-b*t,t)+mq*zerotauint(Cll/2.,b+bl+Cxl,bx+Cxx/2.-b*t,t) +\
        (bq+Cxq-mq)*zerotauint(Cll/2.,b+bl+Cxl-gq,bx+Cxx/2.-b*t,t);
}

double mean_l(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return ml+(bl-ml)*exp(-gl*t);
}
double mean_q(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return mq+(bq-mq)*exp(-gq*t);
}

// ======================================================================================================== //
// COVARIANCE ELEMENTS
// ======================================================================================================== //

double cov_xx(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cll*pow((1-exp(-gl*t)),2)/pow(gl,2)+2*Cxl*(1-exp(-gl*t))/gl+Cxx+ sl2/(2*pow(gl,3))*(2*gl*t-3+4*exp(-gl*t)-pow(exp(-gl*t),2) ) ;
}

double cov_xg(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b, Eigen::VectorXd nm){
	return (bg*bx)/exp(b*t) + Cxg/exp(b*t) + (bg*bl)/(exp(b*t)*gl) + Cgl/(exp(b*t)*gl) - (bg*bl)/(exp((b + gl)*t)*gl) - \
        Cgl/(exp((b + gl)*t)*gl) - (bg*ml)/(exp(b*t)*gl) + (bg*ml)/(exp((b + gl)*t)*gl) + (bg*ml*t)/exp(b*t) + \
        (Cxl*mq + (Cll*mq)/gl)*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t) - \
        (Cll*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t - gl*t,t))/gl + \
        (bx*Clq + bq*Cxl + Cxl*Cxq + Clq*Cxx + (bq*Cll)/gl + (bl*Clq)/gl + (Clq*Cxl)/gl + (Cll*Cxq)/gl - (Clq*ml)/gl - Cxl*mq - \
           (Cll*mq)/gl + Clq*ml*t)*onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t) + \
        (-((bq*Cll)/gl) - (bl*Clq)/gl - (Clq*Cxl)/gl - (Cll*Cxq)/gl + (Clq*ml)/gl + (Cll*mq)/gl)*\
         onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gl*t,t) + \
        (Clq*Cxl + (Cll*Clq)/gl)*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2.- b*t,t) - \
        (Cll*Clq*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gl*t,t))/gl + \
        (bx*mq + Cxx*mq + (bl*mq)/gl + (Cxl*mq)/gl - (ml*mq)/gl + ml*mq*t)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t) + \
        (-((bl*mq)/gl) - (Cxl*mq)/gl + (ml*mq)/gl)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t - gl*t,t) + \
        (bq*bx + Cxq + bx*Cxq + bq*Cxx + Cxq*Cxx + (bl*bq)/gl + Clq/gl + (bq*Cxl)/gl + (bl*Cxq)/gl + (Cxl*Cxq)/gl - (bq*ml)/gl - \
           (Cxq*ml)/gl - bx*mq - Cxx*mq - (bl*mq)/gl - (Cxl*mq)/gl + (ml*mq)/gl + bq*ml*t + Cxq*ml*t - ml*mq*t)*\
         zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t) + \
        (-((bl*bq)/gl) - Clq/gl - (bq*Cxl)/gl - (bl*Cxq)/gl - (Cxl*Cxq)/gl + (bq*ml)/gl + (Cxq*ml)/gl + (bl*mq)/gl + (Cxl*mq)/gl - \
           (ml*mq)/gl)*zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gl*t,t)- nm(1)*nm(0);
}

double cov_xl(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return sl2/(2*pow(gl,2))*pow((1-exp(-gl*t)),2) + Cll*exp(-gl*t)*(1-exp(-gl*t))/gl+Cxl*exp(-gl*t);
}

double cov_xq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Clq*(1-exp(-gl*t))*exp(-gq*t)/gl+Cxq*exp(-gq*t);
}

double cov_gg(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
    return (pow(bg,2) + Cgg)/exp(2*b*t) + \
       2*Cgl*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - 2*b*t,t) + \
       (mq*(2*Clq + gq*mq)*onetauint(Cll/2.,b + bl + 2*Cxl,2*(bx + Cxx - b*t),t))/\
        gq + 2*(bq*Cgl + bg*Clq + Clq*Cxg + Cgl*Cxq - Cgl*mq)*\
        onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - 2*b*t,t) + \
       ((pow(bq,2)*gq + Cqq*gq + 4*bq*Cxq*gq + 4*pow(Cxq,2)*gq - 2*Clq*mq - 2*bq*gq*mq - \
            4*Cxq*gq*mq + gq*pow(mq,2))*\
          onetauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),t))/gq - \
       pow(mq,2)*onetauint(Cll/2.,b + bl + 2*Cxl,2*(bx + Cxx - b*t),2*t,t) - \
       (2*Clq*mq*onetauint(Cll/2.,b + bl + 2*Cxl,2*bx + 2*Cxx - (2*b + gq)*t,2*t,t))/\
        gq - (sq2*onetauint(Cll/2.,b + bl + 2*Cxl - gq,2*bx + 2*Cxx - 2*b*t,t,0))/\
        (2.*gq) + (sq2*onetauint(Cll/2.,b + bl + 2*Cxl - gq,2*bx + 2*Cxx - 2*b*t,2*t,\
           t))/(2.*gq) + (-pow(bq,2) - Cqq - 4*bq*Cxq - 4*pow(Cxq,2) + 2*bq*mq + 4*Cxq*mq - \
          pow(mq,2) + 4*bq*Clq*t + 8*Clq*Cxq*t - 4*Clq*mq*t)*\
        onetauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),2*t,t) + \
       (2*Clq*mq*onetauint(Cll/2.,b + bl + 2*Cxl - gq,2*bx + 2*Cxx - 2*b*t + gq*t,2*t,\
           t))/gq + pow(Clq,2)*treetauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),\
         t) - pow(Clq,2)*treetauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),2*t,\
         t) + 2*Cgl*Clq*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - 2*b*t,t) + \
       (2*bq*Clq + 4*Clq*Cxq - 2*Clq*mq)*\
        twotauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),t) + \
       (-2*bq*Clq - 4*Clq*Cxq + 2*Clq*mq + 2*pow(Clq,2)*t)*\
        twotauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),2*t,t) + \
       (2*bg*mq + 2*Cxg*mq)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - 2*b*t,t) + \
       ((2*bq*mq)/gq + (4*Cxq*mq)/gq - (2*pow(mq,2))/gq)*\
        zerotauint(Cll/2.,b + bl + 2*Cxl,2*(bx + Cxx - b*t),t) + \
       (2*bg*bq + 2*Cgq + 2*bq*Cxg + 2*bg*Cxq + 2*Cxg*Cxq - 2*bg*mq - 2*Cxg*mq)*\
        zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - 2*b*t,t) + \
       ((-2*bq*mq)/gq - (4*Cxq*mq)/gq + (2*pow(mq,2))/gq)*\
        zerotauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),t) + \
       (sq2*zerotauint(Cll/2.,b + bl + 2*Cxl,2*bx + 2*Cxx - 2*b*t,t,0))/(2.*pow(gq,2)) + \
       (sq2*zerotauint(Cll/2.,b + bl + 2*Cxl,2*bx + 2*Cxx - 2*b*t,2*t,t))/\
        (2.*pow(gq,2)) + 2*pow(mq,2)*t*zerotauint(Cll/2.,b + bl + 2*Cxl,2*(bx + Cxx - b*t),\
         2*t,t) + ((-2*bq*mq)/gq - (4*Cxq*mq)/gq + (2*pow(mq,2))/gq)*\
        zerotauint(Cll/2.,b + bl + 2*Cxl,2*bx + 2*Cxx - (2*b + gq)*t,2*t,t) - \
       (sq2*zerotauint(Cll/2.,b + bl + 2*Cxl - gq,2*bx + 2*Cxx - 2*b*t,t,0))/\
        (2.*pow(gq,2)) - (sq2*t*zerotauint(Cll/2.,b + bl + 2*Cxl - gq,\
           2*bx + 2*Cxx - 2*b*t,2*t,t))/gq + \
       (2*pow(bq,2)*t + 2*Cqq*t + 8*bq*Cxq*t + 8*pow(Cxq,2)*t - 4*bq*mq*t - 8*Cxq*mq*t + \
          2*pow(mq,2)*t)*zerotauint(Cll/2.,b + bl + 2*Cxl - gq,2*(bx + Cxx - b*t),2*t,t)\
        + ((2*bq*mq)/gq + (4*Cxq*mq)/gq - (2*pow(mq,2))/gq)*\
        zerotauint(Cll/2.,b + bl + 2*Cxl - gq,2*bx + 2*Cxx - 2*b*t + gq*t,2*t,t) - \
       (sq2*zerotauint(Cll/2.,b + bl + 2*Cxl + gq,2*bx + 2*Cxx - 2*b*t - 2*gq*t,2*t,\
           t))/(2.*pow(gq,2))-pow(nm(1),2);
}

double cov_gl(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
	return (bg*bl)/exp((b + gl)*t) + Cgl/exp((b + gl)*t) + (bg*ml)/exp(b*t) - (bg*ml)/exp((b + gl)*t) + \
        Cll*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t - gl*t,t) + Clq*ml*onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t) + \
        (bq*Cll + bl*Clq + Clq*Cxl + Cll*Cxq - Clq*ml - Cll*mq)*onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gl*t,t) + \
        Cll*Clq*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gl*t,t) + ml*mq*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t) + \
        (bl*mq + Cxl*mq - ml*mq)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t - gl*t,t) + \
        (bq*ml + Cxq*ml - ml*mq)*zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t) + \
        (bl*bq + Clq + bq*Cxl + bl*Cxq + Cxl*Cxq - bq*ml - Cxq*ml - bl*mq - Cxl*mq + ml*mq)*\
         zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gl*t,t) - nm(1)*nm(2);
}

double cov_gq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
	// std::cout << b + bl + Cxl << "\n";
	// std::cout << b + bl + Cxl - gq << "\n";
    // std::cout << b + bl + Cxl + gq << "\n";
    return (bg*bq)/exp((b + gq)*t) + Cgq/exp((b + gq)*t) + (bg*mq)/exp(b*t) - (bg*mq)/exp((b + gq)*t) + \
        Clq*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t - gq*t,t) + Clq*mq*onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t) + \
        (2*bq*Clq + 2*Clq*Cxq - 2*Clq*mq)*onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gq*t,t) + \
        pow(Clq,2)*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gq*t,t) + pow(mq,2)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t) + \
        (bq*mq + Cxq*mq - pow(mq,2))*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t - gq*t,t) + \
        (bq*mq + Cxq*mq - pow(mq,2))*zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t) - \
        (sq2*zerotauint(Cll/2.,b + bl + Cxl - gq,-b*t + bx + Cxx/2. - gq*t,t))/(2.*gq) + \
        (pow(bq,2) + Cqq + 2*bq*Cxq + pow(Cxq,2) - 2*bq*mq - 2*Cxq*mq + pow(mq,2))*zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t - gq*t,t) + \
        (sq2*zerotauint(Cll/2.,b + bl + Cxl + gq,-b*t + bx + Cxx/2. - gq*t,t))/(2.*gq)- nm(1)*nm(3);
}

double cov_ll(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cll*pow(exp(-gl*t),2) + sl2/(2*gl)*(1-pow(exp(-gl*t),2));
}

double cov_lq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return  Clq*exp(-gl*t)*exp(-gq*t);
}

double cov_qq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return sq2/(2*gq)*(1-pow(exp(-gq*t),2)) + Cqq*pow(exp(-gq*t),2);
}

    
void mean_cov_model(MOMAdata &cell, 
                double t, double ml, 
                double gl, double sl2, 
                double mq, double gq, 
                double sq2, double b){
    // Given p(z0)=N(m,C) find p(zt) given no cell division between 0 and t//

    Eigen::VectorXd nm(4);
    Eigen::MatrixXd nC(4, 4);

    double bx=cell.mean(0);
    double bg=cell.mean(1);
    double bl=cell.mean(2); 
    double bq=cell.mean(3);

    double Cxx=cell.cov(0,0);
    double Cxg=cell.cov(0,1);
    double Cxl=cell.cov(0,2);
    double Cxq=cell.cov(0,3);
    double Cgg=cell.cov(1,1);
    double Cgl=cell.cov(1,2);
    double Cgq=cell.cov(1,3);
    double Cll=cell.cov(2,2);
    double Clq=cell.cov(2,3);
    double Cqq=cell.cov(3,3);

    // Mean
    nm(0) = mean_x(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nm(1) = mean_g(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nm(2) = mean_l(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nm(3) = mean_q(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);

    // Cov
    nC(0,1) = nC(1,0) = cov_xg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm);
    nC(0,2) = nC(2,0) = cov_xl(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nC(0,3) = nC(3,0) = cov_xq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);

    nC(1,2) = nC(2,1) = cov_gl(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm);
    nC(1,3) = nC(3,1) = cov_gq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm); 

    nC(2,3) = nC(3,2) = cov_lq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);

    nC(0,0) = cov_xx(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nC(1,1) = cov_gg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm);
    nC(2,2) = cov_ll(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nC(3,3) = cov_qq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    
    cell.mean = nm;
    cell.cov = nC;

    // for(size_t i=0; i<nm.size(); ++i){
    //     if (std::isnan(nm(i))){
    //         std::cout << "NAN-WARNING: nm(i) " << i << "\n";
    //     }
    // }
    // for(size_t i=0; i<nC.rows(); ++i){
    //     for(size_t j=0; j<nC.cols(); ++j){
    //         if (std::isnan(nC(i,j))){
    //             std::cout << "NAN-WARNING: nC(i,j) " << i << " " << j << "\n";
    //         }
    //     }
    // }

}


// ======================================================================================================== //
// CROSS COVARIANCE ELEMENTS
// ======================================================================================================== //

// ============= //
double cross_cxx(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cxx + Cxl*(1-exp(-gl*t))/gl;
}

double cross_cxg(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cxg + Cgl*(1-exp(-gl*t))/gl;
}

double cross_cxl(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cxl + Cll*(1-exp(-gl*t))/gl;
}

double cross_cxq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cxq + Clq*(1-exp(-gl*t))/gl;
}

// ============= //
double cross_cgx(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
    return (bg*bx)/exp(b*t) + Cxg/exp(b*t) + 
        Cxl*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (bx*Clq + bq*Cxl + Cxl*Cxq + Clq*Cxx - Cxl*mq)* \
        onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        Clq*Cxl*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        (bx*mq + Cxx*mq)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (bq*bx + Cxq + bx*Cxq + bq*Cxx + Cxq*Cxx - bx*mq - Cxx*mq)*\
        zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) - nm(1)*bx;
}


double cross_cgg(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
    return pow(bg,2)/exp(b*t) + Cgg/exp(b*t) + \
        Cgl*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (bq*Cgl + bg*Clq + Clq*Cxg + Cgl*Cxq - Cgl*mq)* \
        onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        Cgl*Clq*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        (bg*mq + Cxg*mq)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (bg*bq + Cgq + bq*Cxg + bg*Cxq + Cxg*Cxq - bg*mq - Cxg*mq)* \
        zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) - nm(1)*bg;
}

double cross_cgl(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
    return (bg*bl)/exp(b*t) + Cgl/exp(b*t) + \
        Cll*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (bq*Cll + bl*Clq + Clq*Cxl + Cll*Cxq - Cll*mq)* \
        onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        Cll*Clq*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        (bl*mq + Cxl*mq)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (bl*bq + Clq + bq*Cxl + bl*Cxq + Cxl*Cxq - bl*mq - Cxl*mq)* \
        zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) - nm(1)*bl;
}


double cross_cgq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b,Eigen::VectorXd nm){
    return (bg*bq)/exp(b*t) + Cgq/exp(b*t) + \
        Clq*mq*onetauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (2*bq*Clq + 2*Clq*Cxq - Clq*mq)* \
        onetauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        pow(Clq,2)*twotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) + \
        (bq*mq + Cxq*mq)*zerotauint(Cll/2.,b + bl + Cxl,bx + Cxx/2. - b*t,t,0) + \
        (pow(bq,2) + Cqq + 2*bq*Cxq + pow(Cxq,2) - bq*mq - Cxq*mq)* \
        zerotauint(Cll/2.,b + bl + Cxl - gq,bx + Cxx/2. - b*t,t,0) - nm(1)*bq;
}

// ============= //
double cross_clx(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cxl * exp(-gl*t);
}

double cross_clg(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cgl * exp(-gl*t);
}

double cross_cll(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cll * exp(-gl*t);
}

double cross_clq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Clq * exp(-gl*t);
}

// ============= //
double cross_cqx(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cxq * exp(-gq*t);
}

double cross_cqg(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cgq * exp(-gq*t);
}

double cross_cql(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Clq * exp(-gq*t);
}

double cross_cqq(double t,double bx,double bg,double bl,double bq,double Cxx,double Cxg,double Cxl,double Cxq,double Cgg,double Cgl,double Cgq,double Cll,double Clq,double Cqq,double ml,double gl,double sl2,double mq,double gq,double sq2,double b){
    return Cqq * exp(-gq*t);
}


Eigen::MatrixXd cross_cov_model(MOMAdata &cell, 
                                double t, double ml, 
                                double gl, double sl2, 
                                double mq, double gq, 
                                double sq2, double b){
    Eigen::VectorXd nm(4);
    Eigen::MatrixXd cross_cov(4, 4);

    double bx=cell.mean(0);
    double bg=cell.mean(1);
    double bl=cell.mean(2); 
    double bq=cell.mean(3);

    double Cxx=cell.cov(0,0);
    double Cxg=cell.cov(0,1);
    double Cxl=cell.cov(0,2);
    double Cxq=cell.cov(0,3);
    double Cgg=cell.cov(1,1);
    double Cgl=cell.cov(1,2);
    double Cgq=cell.cov(1,3);
    double Cll=cell.cov(2,2);
    double Clq=cell.cov(2,3);
    double Cqq=cell.cov(3,3);

    // Mean
    nm(0) = mean_x(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nm(1) = mean_g(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nm(2) = mean_l(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);
    nm(3) = mean_q(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b);

    /* 
    <x_n+1, x_n+1>  ...     <x_n+1, x_n>    ...     <q_n+1, q_n> 
    ...
    <q_n, q_n+1>    ...     <q_n, x_n>      ...     <q_n, q_n> 
    */
    cross_cov <<    cross_cxx(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b), // x
                    cross_cxg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b),
                    cross_cxl(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b),
                    cross_cxq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b), 
                    cross_cgx(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm), // g
                    cross_cgg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm),
                    cross_cgl(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm),
                    cross_cgq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b,nm), 
                    cross_clx(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b), // l
                    cross_clg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b),
                    cross_cll(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b),
                    cross_clq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b), 
                    cross_cqx(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b), // q
                    cross_cqg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b),
                    cross_cql(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b),
                    cross_cqq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,b); 
    return cross_cov;
}