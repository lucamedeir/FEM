

#include <lapacke.h>
#include<vector>
#include<cassert>
#include <gtest/gtest.h>
#include <cmath>
#include "builders.hpp"
#include "utils.hpp"
#include "interpol.hpp"
#include <functional>



matrix build_Ce(int n, std::function<double(double)> f) {
    
    assert(n >= 0);

    int size = n+1;
    int width = size;

    double* Ce = new double[size*size];


    for(int row = 0; row < size; row++ ) {
        for(int col = 0; col < size; col++) {

            switch (n)
            {
            case 3: // Quadrático
                Ce[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n3pol(row,x)*n3pol(col,x);},4);
                break;

            case 2: // Quadrático
                Ce[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n2pol(row,x)*n2pol(col,x);},4);
                break;
            
            default: // Linear
                Ce[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n1pol(row,x)*n1pol(col,x);},2);
            }
        }
    }

    matrix result = {Ce, size, size};


    return result;
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

matrix build_Ae(int n, std::function<double(double)> f) {

    assert(n >= 1);

    int size = n+1;
    int width = size;

    double* Ae = new double[size*size];

    for(int row = 0; row < size; row++ ) {
        for(int col = 0; col < size; col++) {

            switch (n)
            {
            case 3: // Cúbico
                Ae[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n3dpol(row,x)*n3dpol(col,x);},4);
                break;

            case 2: // Quadrático
                Ae[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n2dpol(row,x)*n2dpol(col,x);},2);
                break;
            
            default: // Linear
                Ae[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n1dpol(row,x)*n1dpol(col,x);},2);
            }
        }
    }

    matrix result = {Ae, size, size};


    return result;
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


matrix build_Be(int n, std::function<double(double)> f) {

    assert(n >= 1);

    int size = n+1;
    int width = size;

    double* Be = new double[size*size];

    for(int row = 0; row < size; row++ ) {
        for(int col = 0; col < size; col++) {

            switch (n)
            {
            case 3: // Cúbico
                Be[col + row*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n3pol(row,x)*n3dpol(col,x);},4);
                break;

            case 2: // Quadrático
                Be[col + row*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n2pol(row,x)*n2dpol(col,x);},2);
                break;
            
            default: // Linear
                Be[col + row*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n1pol(row,x)*n1dpol(col,x);},2);
            }
        }
    }

    matrix result = {Be, size, size};


    return result;
}

// test Ae
TEST(StiffnessTest, BeTest) {
    auto Be = build_Be(1, [](double x){return 1;});
    std::vector<double> BeTest = {-1.0/2.0, 1.0/2.0, -1.0/2.0, 1.0/2.0}; 

    for(int i = 0; i < Be.nrows*Be.ncols; i++) {
        ASSERT_LT( abs(Be.m[i]-BeTest[i]), 0.00001 ) << "Valores diferentes na matriz Be(1) no índice " << i;
    } 

    delete [] Be.m;
}

vector build_Fe(int n, std::function<double(double)> f) {

    assert(n >= 1);

    int size = n + 1;

    double* Fe = new double[size];

    for(int row = 0; row < size; row++) {
        switch (n)
            {
            case 3: // Cúbico
                Fe[row] = gauss_legendre([&row, f](double x){return f(x)*n3pol(row,x);},3);
                break;

            case 2: // Quadrático
                Fe[row] = gauss_legendre([&row, f](double x){return f(x)*n2pol(row,x);},2);
                break;
            
            default: // Linear
                Fe[row] = gauss_legendre([&row, f](double x){return f(x)*n1pol(row,x);},2);
            }
    }

    vector result = {Fe, size};

    return result;
}

// test Fe
TEST(LoadTest, Fe1Test) {
    auto Fe1 = build_Fe(1,[](double x){return 1;});
    std::vector<double> Fe1Test = {1.0,1.0}; 

    for(int i = 0; i < Fe1.n; i++) {
        ASSERT_LT( abs(Fe1.m[i]- Fe1Test[i]), 0.00001 ) << "Valores diferentes na matriz Fe(1) no índice " << i;
    } 

    delete [] Fe1.m;
}

// test Fe
TEST(LoadTest, Fe3Test) {
    auto Fe2 = build_Fe(3,[](double x){return 1;});
    std::vector<double> Fe2Test = { 0.25, 0.75, 0.75, 0.25 }; 

    for(int i = 0; i < Fe2.n; i++) {
        ASSERT_LT( abs(Fe2.m[i] - Fe2Test[i]), 0.00001 ) << "Valores diferentes na matriz Fe(3) no índice " << i;
    } 

    delete [] Fe2.m;
}

TEST(LoadTest, Fe2Test) {
    auto Fe2 = build_Fe(2,[](double x){return 1;});
    std::vector<double> Fe2Test = {1.0/3.0,4.0/3.0,1.0/3.0}; 

    for(int i = 0; i < Fe2.n; i++) {
        ASSERT_LT( abs(Fe2.m[i] - Fe2Test[i]), 0.00001 ) << "Valores diferentes na matriz Fe(2) no índice " << i;
    } 

    delete [] Fe2.m;
}


std::vector<std::vector<int>> build_LM(int nX, int n) {

    assert(n >= 1);

    int nE = nX - 1;

    std::vector<std::vector<int>> LM;
    LM.reserve(nE);

    
    int index_p = 0;
    for(int i = 0; i < nE; i++) {

        std::vector<int> global_p;
        global_p.reserve(n+1);
        for(int k =0 ; k < n+1; k++) {
            global_p.push_back(index_p + k);
        }
        LM.push_back(global_p);

        index_p += n;
    }

    return LM;
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

matrix build_K(std::vector<double> X, std::vector<std::vector<int>> LM, std::function<double(double)> alpha, std::function<double(double)> beta, std::function<double(double)> gamma, int n) {
    int N = X.size();

    // Número de pontos para K
    // Ne = N-1 (número de elementos)
    // np = n+1 (número de pontos por elementos)
    // Número de pontos = Ne*np -(Ne-1)
    // pois tem que retirar os pontos que conectam os elementos
    int size_K = (N-1)*(n+1) - (N-2);

    double he = X[1]-X[0];

    auto Ae = build_Ae(n,alpha);
    auto Be = build_Ae(n,beta);
    auto Ce = build_Ce(n,gamma);

    double* K = new double[size_K*size_K];
    for(int i = 0; i < size_K*size_K; i++) K[i] = 0;

    int e_index = 0 ; // index do elemento
    int local_i, local_j;
    for(auto e: LM) {

        for(auto global_i: e) {
            for(auto global_j: e) {
                local_i = global_i - n*e_index;
                local_j = global_j - n*e_index;


                K[global_i + global_j*size_K] += 2/he*Ae.m[local_i + local_j*Ae.nrows] + he/2*Ce.m[local_i + local_j*Ce.nrows] + Be.m[local_i + local_j*Ce.nrows];
            }
        }
        e_index++;
    }

    K[0] += 2/he*Ae.m[Ae.nrows*Ae.ncols-1] + he/2*Ce.m[Ce.nrows*Ce.ncols-1] + Be.m[Be.ncols*Be.nrows-1];
    K[size_K*size_K-1] += 2/he*Ae.m[0] + he/2*Ce.m[0] + Be.m[0];


    matrix result = {K, size_K, size_K};

    delete[] Ae.m;
    delete[] Ce.m;

    return result;
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

vector build_F(std::vector<double> X, std::vector<std::vector<int>> LM, std::function<double(double)> f,  int n) {
    int N = X.size();

    // Número de pontos para K
    // Ne = N-1 (número de elementos)
    // np = n+1 (número de pontos por elementos)
    // Número de pontos = Ne*np -(Ne-1)
    // pois tem que retirar os pontos que conectam os elementos
    int sizeF = (N-1)*(n+1) - (N-2);

    double he = X[1]-X[0];

    auto Fe = build_Fe(n,f);

    double* F = new double[sizeF];
    for(int i = 0; i < sizeF; i++) F[i] = 0;

    int e_index = 0 ; // index do elemento
    int local_i;
    for(auto e: LM) {

        for(auto global_i: e) {
            local_i = global_i - n*e_index;


            F[global_i] += he/2*Fe.m[local_i];
        }
        e_index++;
    }

    F[0] += he/2*Fe.m[Fe.n-1];
    F[sizeF-1] += he/2*Fe.m[0];

    vector result = {F, sizeF};

    delete[] Fe.m;

    return result;
}

TEST(BuildingTest, FTest) {
    int N = 2;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,1);

    auto F = build_F(X, LM, [](double x){return 1;}, 1 );

    std::vector<int> FTest = {1,1};

    for(int i = 0; i < F.n; i++) {
        ASSERT_EQ( (int)round(F.m[i]), FTest[i] ) << "Valores diferentes na matriz F(2) no índice " << i;
    }

    N = 3;
    X = linspace(0.0, 1.0, N);
    LM = build_LM(N,1);

    auto F2 = build_F(X, LM, [](double x){return 1;}, 1);

    std::vector<double> FTest3 = {0.5,0.5,0.5};

    for(int i = 0; i < F2.n; i++) {
        ASSERT_LT( abs(F2.m[i]- FTest3[i]),0.00001 ) << "Valores diferentes na matriz F(3) no índice " << i;
    }

    delete[] F.m;
    delete[] F2.m;
}

TEST(BuildingTest, F2Test) {
    int N = 3;
    auto X = linspace(0.0, 1.0, N);
    auto LM = build_LM(N,2);

    auto F = build_F(X, LM,  [](double x){return 1;}, 2);

    std::vector<double> FTest = {1.0/6.0, 1.0/3.0, 1.0/6.0, 1.0/3.0, 1.0/6.0};

    for(int i = 0; i < F.n; i++) {
        ASSERT_LT( abs(F.m[i]-FTest[i]), 0.0001 ) << "Valores diferentes na matriz F(3) no índice " << i;
    }

    delete[] F.m;
}

vector solve(matrix K, vector F) {
    int nrows = K.nrows;
    int ncols = K.ncols;

    double* u = new double[nrows];

    for(int i = 0 ; i < nrows; i++) {
        u[i] = 0.0;
    }

    auto info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', nrows, ncols, 1, K.m, ncols, F.m, 1);

    for(int i = 0 ; i < nrows; i++) {
        u[i] = F.m[i];
    }

    vector result = {u, nrows};

    return result;
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

void set_boundary(matrix K,vector F, double ni, double nf, double ui, double uf) {
    int i = 0;
    int j = 0;

    K.m[0] = 1;
    K.m[K.nrows*K.ncols-1] = 1;

    F.m[0] = ui;
    F.m[F.n-1] = uf;

    /*
    for i in range(n):
        F[i+1] -= ui*K[i+1,0]
        F[N-2-i] -= uf*K[N-2-i,N-1]
        K[i+1,0] = K[0,i+1] = K[N-2-i,N-1] = K[N-1,N-2-i] = 0
    */

   for(int i = 1 ; i < ni+1; i++) {
       F.m[i] -= ui*K.m[i];
       K.m[i] = 0;
       K.m[i*K.nrows] = 0;
   }

   for(i = K.ncols*K.nrows-2 ; i > K.ncols*K.nrows - 1 - (nf + 1) ; i--) {
       F.m[i - (K.nrows-1)*K.ncols ] -= uf*K.m[i];
       K.m[i] = 0;
   }

   i = K.nrows - 1;
   for(int j = K.ncols-2 ; j > K.ncols - 1 - (nf + 1) ; j--) {
       K.m[i + j*K.nrows] = 0;
   }


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