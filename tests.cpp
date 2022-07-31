#include "builders.hpp"
#include "utils.hpp"
#include "interpol.hpp"
#include <gtest/gtest.h>

TEST(Build2DSuite, build_Fe2DTest) {
    std::vector<double> Fe2DTest = {1.0,1.0,1.0,1.0};

    auto Fe2D = build_Fe2D(0,0,0,0,[](double e, double r){return 1.0;});

    for(int i = 0; i < 4; i++) {
        ASSERT_LT(abs(Fe2D.m[i]-Fe2DTest[i]), 0.0001) << "index " << i;
    }

    delete[] Fe2D.m;
}

TEST(Build2DSuite, buildLM2DTestNe3) {
    std::vector<std::vector<int>> LM2DTest = { {1, 2, 5, 6},
                                               {2, 3, 6, 7},
                                               {3, 4, 7, 8},

                                               {5, 6, 9, 10},
                                               {6, 7, 10, 11},
                                               {7, 8, 11, 12},

                                               {9, 10, 13, 14},
                                               {10, 11, 14, 15},
                                               {11, 12, 15, 16}};
    
    

    int nos = 4;
    int nE = 3;

    auto LM2D = build_LM2D(nE);

    ASSERT_EQ(LM2D.size(), LM2DTest.size()) << "Sizes not equal";

    for(int i = 0; i < LM2D.size(); i++) {
        auto e = LM2D[i];
        auto eTest = LM2D[i];

        for(int k = 0; k < nos; k++) {
            ASSERT_EQ(e[k], eTest[k]) << "Error in index " << k;
        }
    }
}

TEST(Build2DSuite, buildLM2DTestNe2) {
    std::vector<std::vector<int>> LM2DTest = { {1, 2, 4, 5},
                                               {2, 3, 5, 6},

                                               {4, 5, 7, 8},
                                               {5, 6, 8, 9}};
    
    

    int nos = 4;
    int nE = 2;

    auto LM2D = build_LM2D(nE);

    ASSERT_EQ(LM2D.size(), LM2DTest.size()) << "Sizes not equal";

    for(int i = 0; i < LM2D.size(); i++) {
        auto e = LM2D[i];
        auto eTest = LM2D[i];

        for(int k = 0; k < nos; k++) {
            ASSERT_EQ(e[k], eTest[k]) << "Error in index " << k;
        }
    }
}

TEST(Build2DSuite, buildCe2DTest) {

    auto Ce = build_Ce2D([](double x, double y){return 1;});
    std::vector<double> CeTest = {4.0, 2.0, 2.0, 1.0, 2.0, 4.0, 1.0, 2.0, 2.0, 1.0, 4.0, 2.0, 1.0, 2.0, 2.0, 4.0};

    for(int i = 0; i < Ce.nrows*Ce.ncols; i++) {
        ASSERT_LT( abs(Ce.m[i]*9.0-CeTest[i]), 0.0001 ) << "Problem at index " << i;
    }

    delete[] Ce.m;
}

TEST(Build2DSuite, buildAe2DTest) {

    auto Ae = build_Ae2D([](double x, double y){return 1;});
    std::vector<double> AeTest = {4.0, -1.0, -1.0, -2.0, -1.0, 4.0, -2.0, -1.0, -1.0, -2.0, 4.0, -1.0, -2.0, -1.0, -1.0, 4.0};

    for(int i = 0; i < Ae.nrows*Ae.ncols; i++) {
        ASSERT_LT( abs(Ae.m[i]*6.0-AeTest[i]), 0.0001 ) << "Problem at index " << i;
    }

    delete[] Ae.m;
}

TEST(Integration2DSuite, GaussLegendre2DTest) {
    ASSERT_LT( abs( gauss_legendre2D(pol2d_a, 3) - 1.0 ), 0.0001 ) << "Integration of Pol2Da";
    ASSERT_LT( abs( gauss_legendre2D([](double x, double y){ return pol2d_a(x,y)*pol2d_a(x,y); }, 5) - 4.0/9.0 ), 0.0001 ) << "Integration of Pol2Da";
}

TEST(Interpolation2DSuite, Pol2Da_Test) {
    ASSERT_LT( abs( pol2d_a(0.0,0.0) - 1.0/4.0 ), 0.0001 ) << " Valor diferente para (x,y) = (0,0)";
    ASSERT_LT( abs( pol2d_a(0.2,0.78) - 0.044 ), 0.0001 ) << " Valor diferente para (x,y) = (0.2,0.78)";
    ASSERT_LT( abs( pol2d_a(-0.76,0.056) - 0.41536 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
}

TEST(Interpolation2DSuite, Pol2Db_Test) {
    ASSERT_LT( abs( pol2d_b(-0.76,0.056) - 0.05664 ), 0.00001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
}

TEST(Interpolation2DSuite, Pol2Dc_Test) {
    ASSERT_LT( abs( pol2d_c(-0.76,0.056) - 0.46464 ), 0.00001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
}

TEST(Interpolation2DSuite, Pol2Dd_Test) {
    ASSERT_LT( abs( pol2d_d(-0.76,0.056) - 0.06336 ), 0.00001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
}

TEST(Interpolation2DSuite, Pol2D_Test) {
    ASSERT_LT( abs( n1pol2D(0,-0.76,0.056) - 0.41536 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
    ASSERT_LT( abs( n1pol2D(1,-0.76,0.056) - 0.05664 ), 0.00001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
    ASSERT_LT( abs( n1pol2D(2,-0.76,0.056) - 0.46464 ), 0.00001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
    ASSERT_LT( abs( n1pol2D(3,-0.76,0.056) - 0.06336 ), 0.00001 ) << " Valor diferente para (x,y) = (-0.76,0.056) ";
}

TEST(Interpolation2DSuite, derivate_Pol2Da_Test) {
    ASSERT_LT( abs( dx_pol2d_a(-0.76,0.056)  + 0.236 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
    ASSERT_LT( abs( dy_pol2d_a(-0.76,0.056) + 0.44 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
}

TEST(Interpolation2DSuite, derivate_Pol2Db_Test) {
    ASSERT_LT( abs( dx_pol2d_b(-0.26,0.156)  - 0.211 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.156)";
    ASSERT_LT( abs( dy_pol2d_b(-0.26,0.056) + 0.185 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.056)";
}

TEST(Interpolation2DSuite, derivate_Pol2Dc_Test) {
    ASSERT_LT( abs( dx_pol2d_c(-0.26,0.156)  + 0.289 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.156)";
    ASSERT_LT( abs( dy_pol2d_c(-0.26,0.056) - 0.315 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.056)";
}

TEST(Interpolation2DSuite, derivate_Pol2Dd_Test) {
    ASSERT_LT( abs( dx_pol2d_d(-0.76,0.056)  - 0.264 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
    ASSERT_LT( abs( dy_pol2d_d(-0.76,0.056) - 0.06 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
}

TEST(Interpolation2DSuite, derivate_dy_Pol2D_Test) {
    ASSERT_LT( abs( n1dy_pol2D(0,-0.76,0.056) + 0.44 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
    ASSERT_LT( abs( n1dy_pol2D(1,-0.26,0.056) + 0.185 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.056)";
    ASSERT_LT( abs( n1dy_pol2D(2,-0.26,0.056) - 0.315 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.056)";
    ASSERT_LT( abs( n1dy_pol2D(3,-0.76,0.056) - 0.06 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
}

TEST(Interpolation2DSuite, derivate_dx_Pol2D_Test) {
    ASSERT_LT( abs( n1dx_pol2D(0,-0.76,0.056)  + 0.236 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
    ASSERT_LT( abs( n1dx_pol2D(1,-0.26,0.156)  - 0.211 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.156)";
    ASSERT_LT( abs( n1dx_pol2D(2,-0.26,0.156)  + 0.289 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.26,0.156)";
    ASSERT_LT( abs( n1dx_pol2D(3,-0.76,0.056)  - 0.264 ), 0.0001 ) << " Valor diferente para (x,y) = (-0.76,0.056)";
}


// test Ce 
TEST(StiffnessTest, Ce1Test) {
    auto Ce1 = build_Ce(1, [](double x){return 1;});
    std::vector<int> Ce1Test = {2,1, 1,2}; // tem que multiplicar por 12 os valores de Ce1.m

    for(int i = 0; i < Ce1.nrows*Ce1.ncols; i++) {
        ASSERT_LT( abs(Ce1.m[i]*3 - Ce1Test[i]), 0.00001 ) << "Valores diferentes na matriz Ce no índice " << i;
    } 

    delete [] Ce1.m;
}

TEST(StiffnessTest, Ce2Test) {
    auto Ce2 = build_Ce(2, [](double x){return 1;});
    std::vector<int> Ce2Test = {4, 2, -1,  2, 16, 2, -1, 2, 4}; // tem que multiplicar por 15 os valores de Ce1.m

    for(int i = 0; i < Ce2.nrows*Ce2.ncols; i++) {
        ASSERT_LT( abs(Ce2.m[i]*15 - Ce2Test[i]), 0.0001  ) << "Valores diferentes na matriz Ce no índice " << i;
    } 

    delete [] Ce2.m;
}


// test Ae
TEST(StiffnessTest, Ae1Test) {
    auto Ae1 = build_Ae(1, [](double x){return 1;});
    std::vector<int> Ae1Test = {1,-1, -1,1}; 

    for(int i = 0; i < Ae1.nrows*Ae1.ncols; i++) {
        ASSERT_EQ( (int)round(Ae1.m[i]*2), Ae1Test[i] ) << "Valores diferentes na matriz Ae no índice " << i;
    } 

    delete [] Ae1.m;
}

// test Ae
TEST(StiffnessTest, Ae2Test) {
    auto Ae2 = build_Ae(2, [](double x){return 1;});
    std::vector<int> Ae2Test = {7, -8, 1, -8, 16, -8, 1, -8, 7}; 

    for(int i = 0; i < Ae2.nrows*Ae2.ncols; i++) {
        ASSERT_EQ( (int)round(Ae2.m[i]*6), Ae2Test[i] ) << "Valores diferentes na matriz Ae no índice " << i;
    } 

    delete [] Ae2.m;
}


// test Ae
TEST(StiffnessTest, BeTest) {
    auto Be = build_Be(1, [](double x){return 1;});
    std::vector<double> BeTest = {-1.0/2.0, -1.0/2.0, 1.0/2.0, 1.0/2.0}; 

    for(int i = 0; i < Be.nrows*Be.ncols; i++) {
        ASSERT_LT( abs(Be.m[i]-BeTest[i]), 0.00001 ) << "Valores diferentes na matriz Be(1) no índice " << i;
    } 

    delete [] Be.m;
}


// test Fe
TEST(LoadTest, Fe1Test) {
    auto Fe1 = build_Fe(1,0,1,[](double x){return 1;});
    std::vector<double> Fe1Test = {1.0,1.0}; 

    for(int i = 0; i < Fe1.n; i++) {
        ASSERT_LT( abs(Fe1.m[i]- Fe1Test[i]), 0.00001 ) << "Valores diferentes na matriz Fe(1) no índice " << i;
    } 

    delete [] Fe1.m;
}

// test Fe
TEST(LoadTest, Fe3Test) {
    auto Fe2 = build_Fe(3,0,1,[](double x){return 1;});
    std::vector<double> Fe2Test = { 0.25, 0.75, 0.75, 0.25 }; 

    for(int i = 0; i < Fe2.n; i++) {
        ASSERT_LT( abs(Fe2.m[i] - Fe2Test[i]), 0.00001 ) << "Valores diferentes na matriz Fe(3) no índice " << i;
    } 

    delete [] Fe2.m;
}

TEST(LoadTest, Fe2Test) {
    auto Fe2 = build_Fe(2,0,1,[](double x){return 1;});
    std::vector<double> Fe2Test = {1.0/3.0,4.0/3.0,1.0/3.0}; 

    for(int i = 0; i < Fe2.n; i++) {
        ASSERT_LT( abs(Fe2.m[i] - Fe2Test[i]), 0.00001 ) << "Valores diferentes na matriz Fe(2) no índice " << i;
    } 

    delete [] Fe2.m;
}


// test LM
TEST(BuildingTest, LMTest) {
    auto LM = build_LM(7,1);

    int i = 0;
    int element_i = 1;
    for(auto e: LM) {
        ASSERT_EQ( e[0], i) << "Ponto 1 do Valor do elemento " << element_i;
        ASSERT_EQ( e[1], i+1) << "Ponto 2 do Valor do elemento " << element_i;
        element_i++;
        i+=1;
    }

}

// test LM
TEST(BuildingTest, LM2Test) {
    auto LM = build_LM(7,2);

    int i = 0;
    int element_i = 1;
    for(auto e: LM) {
        ASSERT_EQ( e[0], i) << "Ponto 1 do Valor do elemento " << element_i;
        ASSERT_EQ( e[1], i+1) << "Ponto 2 do Valor do elemento " << element_i;
        ASSERT_EQ( e[2], i+2) << "Ponto 2 do Valor do elemento " << element_i;
        element_i++;
        i+=2;
    }

}



// test K
TEST(BuildingTest, KTest) {

    int N = 2;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);

    auto alpha = [](double){return 1;};
    auto beta = [](double){return 0;};
    auto gamma = [](double){return 0;};

    auto K = build_K(X, LM, alpha, beta, gamma, 1);

    std::vector<int> KTest = {2,-1,-1,2};

    for(int i = 0; i < K.nrows*K.ncols; i++) {
        ASSERT_EQ( (int)round(K.m[i]), KTest[i] ) << "Valores diferentes na matriz K(2) no índice " << i;
    }

    N = 3;
    X = linspace(0.0, 1.0, N);
    LM = build_LM(N,1);

    auto K2 = build_K(X, LM, alpha,beta , gamma,  1);

    std::vector<int> KTest2 = {4,-2,0,-2,4,-2,0,-2,4};

    for(int i = 0; i < K2.nrows*K2.ncols; i++) {
        ASSERT_EQ( (int)round(K2.m[i]), KTest2[i] ) << "Valores diferentes na matriz K(3) no índice " << i;
    }

    delete[] K.m;
    delete[] K2.m;
}

TEST(BuildingTest, K2Test) {

    int N = 3;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,2);

    auto alpha = [](double){return 1;};
    auto beta = [](double){return 0;};
    auto gamma = [](double){return 0;};

    auto K = build_K(X, LM, alpha, beta, gamma, 2);

    std::vector<double> KTest = {9.33333333, -5.33333333,  0.66666667,  0.        ,  0.,
                                -5.33333333, 10.66666667, -5.33333333,  0.        ,  0. ,
                                 0.66666667, -5.33333333,  9.33333333, -5.33333333,  0.66666667,
                                 0.        ,  0.        , -5.33333333, 10.66666667, -5.33333333,
                                 0.        ,  0.        ,  0.66666667, -5.33333333,  9.33333333};

    for(int i = 0; i < K.nrows*K.ncols; i++) {
        ASSERT_LT( abs(K.m[i]-KTest[i]), 0.00001 ) << "Valores diferentes na matriz K(3) no índice " << i;
    }

    delete[] K.m;
}


TEST(BuildingTest, FTest) {
    int N = 2;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);

    auto F = build_F(X, LM, [](double x){return 1;}, 1 );

    std::vector<double> FTest = {1.0,1.0};

    for(int i = 0; i < F.n; i++) {
        ASSERT_LT( abs(F.m[i]-FTest[i]), 0.0001 ) << "F["<<i<<"] = " << F.m[i] << " || FTest["<<i<<"] = "<<FTest[i];
    }

    delete[] F.m;
}

TEST(BuildingTest, FTestN2) {
    int N = 3;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);

    auto F = build_F(X, LM, [](double x){return 1;}, 1 );

    std::vector<double> FTest3 = {0.5,0.5,0.5};

    for(int i = 0; i < F.n; i++) {
        ASSERT_LT( abs(F.m[i]-FTest3[i]), 0.0001 )  << "F["<<i<<"] = " << F.m[i] << " || FTest["<<i<<"] = "<<FTest3[i];
    }

    delete[] F.m;
}

TEST(BuildingTest, F2Test) {
    int N = 3;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,2);

    auto F = build_F(X, LM,  [](double x){return 1;}, 2);

    std::vector<double> FTest = {1.0/6.0, 1.0/3.0, 1.0/6.0, 1.0/3.0, 1.0/6.0};

    for(int i = 0; i < F.n; i++) {
        ASSERT_LT( abs(F.m[i]-FTest[i]), 0.0001 ) << "F["<<i<<"] = " << F.m[i] << " || FTest["<<i<<"] = "<<FTest[i];
    }

    delete[] F.m;
}


TEST(solverTest,solver) {
    int N = 2;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);

    auto alpha = [](double){return 1;};
    auto beta = [](double){return 0;};
    auto gamma = [](double){return 0;};
    auto f = [](double x){return 1;};

    auto K = build_K(X, LM, alpha, beta, gamma,1);
    auto F = build_F(X, LM, f, 1);

    auto u = solve(K,F);

    std::vector<double> uTest = {1.0,1.0};

    for(int i = 0; i < u.n; i++) {
        ASSERT_LT( abs(u.m[i]-uTest[i]),0.0001) << "Test de u(2) na linha i" << i;
    }


    N = 3;
    X = linspace(0.0, 1.0, N);
    LM = build_LM(N,1);
    auto K2 = build_K(X, LM, alpha,beta,gamma,1);
    auto F2 = build_F(X, LM, f, 1);

    auto u2 = solve(K2,F2);

    std::vector<double> uTest2 = {0.375, 0.5  , 0.375};

    for(int i = 0; i < u2.n; i++) {
        ASSERT_LT( abs( u2.m[i] - uTest2[i] ) , 0.0001) << "Test de u(3) na linha i" << i;
    }

    delete[] u.m;
    delete[] K.m;
    delete[] F.m;
    delete[] u2.m;
    delete[] K2.m;
    delete[] F2.m;
}


TEST(solverTest,boundary) {
    int N = 3;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);


    auto alpha = [](double){return 1;};
    auto beta = [](double){return 0;};
    auto gamma = [](double){return 0;};
    auto f = [](double x){return 1;};

    auto K = build_K(X, LM, alpha,beta, gamma,1);
    auto F = build_F(X, LM, f,1);

    set_boundary(K,F,1,1,0,0);

    std::vector<double> FTest = {0,0.5,0};
    std::vector<double> KTest = {1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 1.0};

    for(int i = 0; i < F.n; i++) {
        ASSERT_LT( abs(F.m[i] - FTest[i]), 0.00001 ) << "Test de F(3) na linha i" << i;
    }

    for(int i = 0; i < K.nrows*K.ncols; i++) {
        ASSERT_LT( abs(K.m[i] - KTest[i]), 0.00001 ) << "Test de K(3) na linha i" << i;
    }

    delete[] K.m;
    delete[] F.m;
}

TEST(solverTest,boundary_nonHomogeneous) {
    int N = 5;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);


    auto alpha = [](double){return 1;};
    auto beta = [](double){return 0;};
    auto gamma = [](double){return 0;};
    auto f = [](double x){return 1;};
    
    auto K = build_K(X, LM, alpha,beta,gamma,1);
    auto F = build_F(X, LM, f, 1);

    set_boundary(K,F,1,1,1.2,1.6); // ui=1.2 e uf=1.6

    std::vector<double> FTest = {1.2 , 5.05, 0.25, 6.65, 1.6  };

    for(int i = 0; i < F.n; i++) {
        ASSERT_LT( abs(F.m[i] - FTest[i]), 0.00001 ) << "Test de F(3) na linha i" << i;
    }

    delete[] K.m;
    delete[] F.m;
}


// test gauss_legendre 
TEST(IntegralTest, gauss_legendreTest) {

    auto result = gauss_legendre([](double x){return x*x;},2);

    ASSERT_EQ((int)round(result*3), 2) << "Teste da integral de x*x";

    result = gauss_legendre([](double x){return x*x*x*x;},4);

    ASSERT_EQ((int)round(result*5), 2) << "Teste da integral de x*x*x*x";

   auto f = [](double x){return x*x -x + 1.0/4.0;};
   result = gauss_legendre(f,2);

    ASSERT_EQ((int)round(result*6), 7) << "Teste da integral de x*x - x + 1/4";

    

}


TEST(BuildingTest, linspaceTest) {

    auto X = linspace(0.0, 1.0, 5);

    std::vector<double> XTest = {0.0, 0.25, 0.5, 0.75, 1.0};

    int i = 0;
    for(auto v: X) {
        ASSERT_EQ(v,XTest[i++]) << "position i in vector X:" << i;
    }
}