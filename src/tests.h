
// Naive implementations of critical integrals
double zerotauint2(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 exp[a*s**2+b*s+c]ds//
    double x = (exp(-pow(b,2.)/(4.*a) + c)*sqrt(M_PI)*(-Faddeeva::erfi((b + 2.*a*t0)/(2.*sqrt(a))) + Faddeeva::erfi((b + 2.*a*t1)/(2.*sqrt(a)))))/(2.*sqrt(a));
    if (std::isnan(x)){
        std::cout   << "INF-WARNING: zerotauint (a,b,c,t0,t1) " 
                    << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "  ";
        std::cout << "exp, erfi t0, x, erfi t1, x, "
                    << exp(-pow(b,2)/(4.*a) + c) << " " 
                    << -Faddeeva::erfi((b + 2*a*t0)/(2.*sqrt(a))) << " " << (b + 2*a*t0)/(2.*sqrt(a)) << " " 
                    <<  Faddeeva::erfi((b + 2*a*t1)/(2.*sqrt(a))) << " " << (b + 2*a*t1)/(2.*sqrt(a)) << " " << "\n";
    }
    return x;
}
double onetauint2(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 s*exp[a*s**2+b*s+c]ds//
    double x= (exp(-pow(b,2)/(4.*a) + c)*(-2*sqrt(a)*exp(pow(b,2)/(4.*a))*(exp(t0*(b + a*t0)) - exp(t1*(b + a*t1))) +\
           b*sqrt(M_PI)*Faddeeva::erfi((b + 2*a*t0)/(2.*sqrt(a))) - b*sqrt(M_PI)*Faddeeva::erfi((b + 2*a*t1)/(2.*sqrt(a)))))/(4.*pow(a,1.5));
    if (std::isnan(x)){
        std::cout   << "INF-WARNING: onetauint (a,b,c,t0,t1) " 
                    << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "\n";
    }
    return x;
}

double twotauint2(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 s**2*exp[a*s**2+b*s+c]ds//
    double x = (exp(-pow(b,2)/(4.*a) + c)*(-2*sqrt(a)*exp(pow(b,2)/(4.*a))*\
           (-(b*exp(t0*(b + a*t0))) + b*exp(t1*(b + a*t1)) + 2*a*exp(t0*(b + a*t0))*t0 - 2*a*exp(t1*(b + a*t1))*t1) +\
           (2*a - pow(b,2))*sqrt(M_PI)*Faddeeva::erfi((b + 2*a*t0)/(2.*sqrt(a))) + (-2*a + pow(b,2))*sqrt(M_PI)*Faddeeva::erfi((b + 2*a*t1)/(2.*sqrt(a)))))/(8.*pow(a,2.5));
    if (std::isnan(x)){
        std::cout   << "INF-WARNING: twotauint (a,b,c,t0,t1) " 
                    << a << " "<< b << " "<< c << " "<<  t0 <<  " "<<  t1 << "\n";
    }
    return x;
}



double treetauint2(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 s**3*exp[a*s**2+b*s+c]ds//
    double x = (exp(-pow(b,2)/(4.*a) + c)*(-2*sqrt(a)*exp(pow(b,2)/(4.*a))*\
           (pow(b,2)*(exp(t0*(b + a*t0)) - exp(t1*(b + a*t1))) - 2*a*exp(t0*(b + a*t0))*(2 + b*t0) + 2*a*exp(t1*(b + a*t1))*(2 + b*t1) +\
            4*pow(a,2)*(exp(t0*(b + a*t0))*pow(t0,2) - exp(t1*(b + a*t1))*pow(t1,2))) + b*(-6*a + pow(b,2))*sqrt(M_PI)*Faddeeva::erfi((b + 2*a*t0)/(2.*sqrt(a))) -\
           b*(-6*a + pow(b,2))*sqrt(M_PI)*Faddeeva::erfi((b + 2*a*t1)/(2.*sqrt(a)))))/(16.*pow(a,3.5));
    if (std::isnan(x)){
        std::cout   << "INF-WARNING: treetauint (a,b,c,t0,t1) " 
                    << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "\n";
    }
    return x;
}

double zerotauint_opt(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 exp[a*s**2+b*s+c]ds rewritten with Dawson functions 
    double x = 2. *(- exp(a*pow(t0,2) + b*t0 +c) * Faddeeva::Dawson((b + 2.*a*t0)/(2.*sqrt(a))) \
                    + exp(a*pow(t1,2) + b*t1 +c) * Faddeeva::Dawson((b + 2.*a*t1)/(2.*sqrt(a))) );

    // if (std::isnan(x)){
    //     std::cout   << "INF-WARNING: zerotauint (a,b,c,t0,t1) " << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "  ";
    // }
    return x /(2.*sqrt(a));
}

double onetauint_opt(double a, double b, double c, double t1, double t0=0){
    //int_t0^t1 s*exp[a*s**2+b*s+c]ds rewritten with Dawson functions 
    double x= -2.*sqrt(a) * (exp(c + t0*(b + a*t0)) - exp(c + t1*(b + a*t1))) +\
           b*2.* (exp(-pow(b,2.)/(4.*a) + c + pow((b + 2.*a*t0)/(2.*sqrt(a)), 2.)) * Faddeeva::Dawson((b + 2.*a*t0)/(2.*sqrt(a))) \
           -exp(-pow(b,2.)/(4.*a) + c + pow((b + 2.*a*t1)/(2.*sqrt(a)), 2.)) * Faddeeva::Dawson((b + 2.*a*t1)/(2.*sqrt(a))));

    // if (std::isnan(x)){
    //     std::cout   << "INF-WARNING: onetauint (a,b,c,t0,t1) " << a << " "<< b << " "<< c << " "<<  t0 << " "<<  t1 << "\n";
    // }    
    return x/(4.*pow(a,1.5));
}

double twotauint_opt(double a, double b, double c, double t_1, double t_0=0){
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

void test_mean_cov_model(){
    double a = 0.0111;
    double b = 0.022;
    double c = 0.01;
    double t1 = 0.2;
    // double t0 = 2;
    double t0 = 0.7;


    std::cout << "---------- ***tauint -----------"<< "\n";
    std::cout << "zerotauint     " << zerotauint(a, b, c, t1, t0) << "\n";
    std::cout << "zerotauint2    " << zerotauint2(a, b, c, t1, t0) << "\n";
    std::cout << "zerotauint_opt " << zerotauint_opt(a, b, c, t1, t0) << "\n\n";

    std::cout << "onetauint     " << onetauint(a, b, c, t1, t0) << "\n";
    std::cout << "onetauint2    " << onetauint2(a, b, c, t1, t0) << "\n";
    std::cout << "onetauint_opt " << onetauint_opt(a, b, c, t1, t0) << "\n\n";

    std::cout << "twotauint     " << twotauint(a, b, c, t1, t0) << "\n";
    std::cout << "twotauint2    " << twotauint2(a, b, c, t1, t0) << "\n";
    std::cout << "twotauint_opt " << twotauint_opt(a, b, c, t1, t0) << "\n\n";
    
    std::cout << "treetauint     "  << treetauint(a, b, c, t1, t0) << "\n";
    std::cout << "treetauint2    "  << treetauint2(a, b, c, t1, t0) << "\n";
    std::cout << "treetauint_opt "  << treetauint_opt(a, b, c, t1, t0) << "\n\n";


    double 	t	 = 	0.1	;
    double 	bx 	 = 	0.2	;
    double 	bg 	 = 	0.3	;
    double 	bl 	 = 	0.4	;
    double 	bq 	 = 	0.5	;
    double 	Cxx	 = 	0.6	;
    double 	Cxg	 = 	0.7	;
    double 	Cxl	 = 	0.8	;
    double 	Cxq	 = 	0.9	;
    double 	Cgg	 = 	1	;
    double 	Cgl	 = 	1.1	;
    double 	Cgq	 = 	1.2	;
    double 	Cll	 = 	1.3	;
    double 	Clq	 = 	1.4	;
    double 	Cqq	 = 	1.5	;
    double 	ml 	 = 	1.6	;
    double 	gl 	 = 	1.7	;
    double 	sl2	 = 	1.8	;
    double 	mq 	 = 	1.9	;
    double 	gq 	 = 	2	;
    double 	sq2	 = 	2.1	;
    double 	beta = 	2.2	;

    Eigen::VectorXd nm(4); 

    nm(0) = 1;
    nm(1) = 2;
    nm(2) = 3;
    nm(3) = 4;

    std::cout << "---------- mean-cov terms -----------"<< "\n";
    std::cout << mean_x(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta)<< "\n" ;
    std::cout << mean_g(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << mean_l(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << mean_q(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << cov_xx(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << cov_xg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta,nm) << "\n" ;  
    std::cout << cov_xl(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << cov_xq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << cov_gg(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta,nm) << "\n" ;  
    std::cout << cov_gl(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta,nm) << "\n" ;  
    std::cout << cov_gq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta,nm) << "\n" ;  
    std::cout << cov_ll(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << cov_lq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  
    std::cout << cov_qq(t,bx,bg,bl,bq,Cxx,Cxg,Cxl,Cxq,Cgg,Cgl,Cgq,Cll,Clq,Cqq,ml,gl,sl2,mq,gq,sq2,beta) << "\n" ;  


    /******************************************/
    Eigen::VectorXd xg(2);
    xg(0) = 20;
    xg(1) = 10;

    MOMAdata cell;
    cell.cov(0,0) = 1;
    cell.cov(1,1) = 2;
    cell.cov(2,2) = 3;
    cell.cov(3,3) = 4;

    cell.cov(1,0) = 2;
    cell.cov(0,1) = 2;
    cell.cov(3,1) = 3;
    cell.cov(1,3) = 4;

    cell.mean = nm;

    xg(0) = xg(0) - cell.mean(0);
    xg(1) = xg(1) - cell.mean(1);

    Eigen::MatrixXd D(2,2);
    D <<  5, 0, 0,  5;

    Eigen::Matrix2d S;
    Eigen::Matrix2d Si;

    S = cell.cov.block(0,0,2,2) + D;
    Si = S.inverse();

    std::cout << "---------- MEAN COV before -----------"<< "\n";
    std::cout << cell << "\n";
    std::cout << xg << "\n";

    mean_cov_model(cell, 1 , 1, 
                        2, 3, 4, 
                        5, 6, 7);

    std::cout << "---------- MEAN COV after mean_cov_model -----------"<< "\n";
    std::cout << cell << "\n";

}

// void test_division(){
//     MOMAdata cell;
//     MOMAdata daugther_cell;

//     daugther_cell.parent = &cell;
    
//     Eigen::VectorXd nm(4); 

//     nm(0) = 1;
//     nm(1) = 2;
//     nm(2) = 3;
//     nm(3) = 4;

//     cell.mean = nm;

//     cell.cov(0,0) = 1;
//     cell.cov(1,1) = 2;
//     cell.cov(2,2) = 3;
//     cell.cov(3,3) = 4;

//     cell.cov(1,0) = 2;
//     cell.cov(0,1) = 2;
//     cell.cov(3,1) = 3;
//     cell.cov(1,3) = 3;

//     mean_cov_after_division(daugther_cell, 0.5, 0.5);

//     std::cout << "---------- MEAN COV after division -----------"<< "\n";
//     std::cout << daugther_cell << "\n";
// }


// void test_likelihood(){
//     // Y,m,C

//         MOMAdata cell;
//         cell.cov << 4.25476409e-02,  4.81488709e+01, -6.17116203e-05, -1.25892662e-01, 
//                     4.81488709e+01,  1.67680116e+06,  2.59861605e-01, 7.45274531e+02,
//                     -6.17116203e-05,  2.59861605e-01,  8.48575294e-07, 8.44383560e-05,
//                     -1.25892662e-01,  7.45274531e+02,  8.44383560e-05, 1.63738212e+00;

//         cell.mean <<    6.93147181e-01,
//                         6.03801845e+03,
//                         1.00811380e-02,
//                         9.56031050e+00;
//         cell.log_length.resize(3);
//         cell.log_length << 0.6621376048238568, 0.8057995840040671, 1.016156660637409;

//         cell.fp.resize(3);
//         cell.fp << 6031.236638936349, 6179.754023612084 , 6351.815340631341;

//         cell.time.resize(3);
//         cell.time << 0, 15, 30;

//         std::cout << "---------- LIKELIHOOD -----------"<< "\n";

//         std::cout << "time:\n" << cell.time<< "\nlog_length:\n" << cell.log_length << "\nfp:\n" << cell.fp<< "\n";

//         std::cout << cell << "\n";

//         std::vector<double> params_vec = {0.01,
//                                             0.01,
//                                             1e-05,
//                                             10,
//                                             0.01,
//                                             0.1,
//                                             0.001,
//                                             0.001,
//                                             5000.0,
//                                             0.001,
//                                             5000.0};
//         double tl = 0;
//         sc_likelihood(params_vec, cell, tl);

//         std::cout << tl;
// }


// void test_prediction(){
//     // Y,m,C

//         MOMAdata cell;
//         cell.cov << 4.25476409e-02,  4.81488709e+01, -6.17116203e-05, -1.25892662e-01, 
//                     4.81488709e+01,  1.67680116e+06,  2.59861605e-01, 7.45274531e+02,
//                     -6.17116203e-05,  2.59861605e-01,  8.48575294e-07, 8.44383560e-05,
//                     -1.25892662e-01,  7.45274531e+02,  8.44383560e-05, 1.63738212e+00;

//         cell.mean <<    6.93147181e-01,
//                         6.03801845e+03,
//                         1.00811380e-02,
//                         9.56031050e+00;
//         cell.log_length.resize(3);
//         cell.log_length << 0.6621376048238568, 0.8057995840040671, 1.016156660637409;

//         cell.fp.resize(3);
//         cell.fp << 6031.236638936349, 6179.754023612084 , 6351.815340631341;

//         cell.time.resize(3);
//         cell.time << 0, 15, 30;

//         std::cout << "---------- PREDICTIONS -----------"<< "\n";

//         std::cout << "time:\n" << cell.time<< "\nlog_length:\n" << cell.log_length << "\nfp:\n" << cell.fp<< "\n";

//         std::cout << cell << "\n";

//         std::vector<double> params_vec = {0.01,
//                                             0.01,
//                                             1e-05,
//                                             10,
//                                             0.01,
//                                             0.1,
//                                             0.001,
//                                             0.001,
//                                             5000.0,
//                                             0.001,
//                                             5000.0};
//         sc_prediction_forward(params_vec, cell);
//         for (size_t i =0; i<cell.mean_forward.size();++i){
//           std::cout << cell.mean_forward[i]; 

//         }
//         std::cout << cell;
// }

// void run_likelihood(CSVconfig config, Parameter_set params, std::string infile){

//     std::cout << "-> Reading" << "\n";
//     std::vector<MOMAdata> cells =  read_data(infile, config);

//     /* genealogy */
//     build_cell_genealogy(cells);

//     Eigen::VectorXd mean(4);
//     mean << 6.93147181e-01,
//             6.02556189e+03,
//             1.03989065e-02,
//             1.02454986e+01;

//     Eigen::MatrixXd cov(4,4);
//     cov <<  1.22158419e-04, -3.38642002e-01,  3.42444314e-06, -4.90827026e-04,
//             -3.38642002e-01,  1.25286734e+05, -3.88680250e-01,  1.42591667e+02,
//             3.42444314e-06, -3.88680250e-01,  4.47368172e-06,  5.05127089e-05,
//             -4.90827026e-04,  1.42591667e+02,  5.05127089e-05,  2.38674307e+00;

//     init_cells_f(cells, mean, cov);

//     double tl = 0;
//     pvector(params.get_init());
//     sc_likelihood(params.get_init(), cells[0], tl);
//     std::cout << "tl: " << tl << "\n";
// }


// void init_cells_f(std::vector<MOMAdata> &cells, Eigen::VectorXd mean, Eigen::MatrixXd cov){
//     /* 
//     * Inititalizes the mean vector and the covariance matrix of the root cells with
//     * pre-defined values
//     */
//     for(size_t i=0; i<cells.size(); ++i){
//         cells[i].mean_init_forward = Eigen::VectorXd::Zero(4);
//         cells[i].cov_init_forward = Eigen::MatrixXd::Zero(4, 4);
//     }

//     std::vector<MOMAdata *> roots = get_roots(cells);
//     for(size_t i=0; i<roots.size(); ++i){
//         roots[i]->mean_init_forward = mean;
//         roots[i]->cov_init_forward = cov;
//     }
// }
