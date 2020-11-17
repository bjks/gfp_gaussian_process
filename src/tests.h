
void test_cases(){
    double a = 0.1;
    double b = 0.2;
    double c = 0.3;
    double t1 = 1;
    double t0 = 2;

    std::cout << zerotauint(a, b, c, t1, t0) << "\n";
    std::cout << onetauint(a, b, c, t1, t0) << "\n";
    std::cout << twotauint(a, b, c, t1, t0) << "\n";
    std::cout << treetauint(a, b, c, t1, t0) << "\n";

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


    std::cout << cell.cov << "\n";
    std::cout << cell.mean << "\n";
    std::cout << xg << "\n";


    mean_cov_model(cell, 1 , 1, 
                        2, 3, 4, 
                        5, 6, 7);
    std::cout << cell.cov << "\n";
    std::cout << cell.mean << "\n";               

}

void test_likelihood(){
    // Y,m,C
/*
(array([[0.6621376048238568, 6031.236638936349, 0.0, '00.0'],
        [0.8057995840040671, 6179.754023612084, 15.0, '00.0'],
        [1.016156660637409, 6351.815340631341, 30.0, '00.0']], dtype=object),
 array([[6.93147181e-01],
        [6.03801845e+03],
        [1.00811380e-02],
        [9.56031050e+00]]),
 array([[ 4.25476409e-02,  4.81488709e+01, -6.17116203e-05,
         -1.25892662e-01],
        [ 4.81488709e+01,  1.67680116e+06,  2.59861605e-01,
          7.45274531e+02],
        [-6.17116203e-05,  2.59861605e-01,  8.48575294e-07,
          8.44383560e-05],
        [-1.25892662e-01,  7.45274531e+02,  8.44383560e-05,
          1.63738212e+00]]),

        ml,gl,sl2,mq,gq,sq2,b,sx2,sg2,sdx2,sdg2
        0.01,
        0.01,
        1e-07,
        10,
        0.01,
        0.1,
        0.001,
        0.001,
        5000.0,
        0.001,
        5000.0)
        */
        MOMAdata cell;
        cell.cov << 4.25476409e-02,  4.81488709e+01, -6.17116203e-05, -1.25892662e-01, 
                    4.81488709e+01,  1.67680116e+06,  2.59861605e-01, 7.45274531e+02,
                    -6.17116203e-05,  2.59861605e-01,  8.48575294e-07, 8.44383560e-05,
                    -1.25892662e-01,  7.45274531e+02,  8.44383560e-05, 1.63738212e+00;

        cell.mean <<    6.93147181e-01,
                        6.03801845e+03,
                        1.00811380e-02,
                        9.56031050e+00;
        cell.log_length.resize(3);
        cell.log_length << 0.6621376048238568, 0.8057995840040671, 1.016156660637409;

        cell.fp.resize(3);
        cell.fp << 6031.236638936349, 6179.754023612084 , 6351.815340631341;

        cell.time.resize(3);
        cell.time << 0, 15, 30;


        std::cout << "time " << cell.time<< "\nlog_length " << cell.log_length << "\nfp" << cell.fp<< "\n";

        std::cout << cell.mean<< "\n" << cell.cov<< "\n";

        std::vector<double> params_vec = {0.01,
                                            0.01,
                                            1e-05,
                                            10,
                                            0.01,
                                            0.1,
                                            0.001,
                                            0.001,
                                            5000.0,
                                            0.001,
                                            5000.0};
        double tl = 0;
        sc_likelihood(params_vec, cell, tl);

        std::cout << tl;
}
