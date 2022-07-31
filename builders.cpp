

#include <lapacke.h>
#include<vector>
#include<cassert>
#include <cmath>
#include "builders.hpp"
#include "utils.hpp"
#include "interpol.hpp"
#include <functional>


vector build_Fe2D(double ax, double bx, double ay, double by, std::function<double(double,double)> f) {

    int size = 4;

    double* Fe = new double[size];

    for(int row = 0; row < size; row++) {
        Fe[row] = gauss_legendre2D([&row, &ax, &bx, &ay, &by, f](double e, double r){return f( (bx-ax)/2*(e+1)+ax, (by-ay)/2*(r+1)+ay )*n1pol2D(row,e,r);},5);
    }

    vector result = {Fe, size};

    return result;
}

matrix build_K2D(std::vector<double> X,
               std::vector<std::vector<int>> LM, 
               std::function<double(double,double)> alpha, 
               std::function<double(double,double)> gamma) {
    int N = X.size();

    int size_K = N;

    double he = X[1]-X[0];

    auto Ae = build_Ae2D(alpha);
    auto Ce = build_Ce2D(gamma);

    double* K = new double[size_K*size_K];
    for(int i = 0; i < size_K*size_K; i++) K[i] = 0;

    int e_index = 0 ; // index do elemento
    int local_i, local_j;
    for(auto e: LM) {

        for(int local_i = 0; local_i < e.size(); local_i ++) {
            for(int local_j = 0; local_j < e.size(); local_j ++) {
                auto global_i = e[local_i]-1;
                auto global_j = e[local_j]-1;

                K[global_i + global_j*size_K] += 4/(he*he)*Ae.m[local_i + local_j*Ae.nrows] + he*he/4*Ce.m[local_i + local_j*Ce.nrows];
            }
        }
    }

    matrix result = {K, size_K, size_K};

    delete[] Ae.m;
    delete[] Ce.m;

    return result;
}

std::vector<std::vector<int>> build_LM2D(int nE) {


    std::vector<std::vector<int>> LM;

    const int nEX = nE;
    const int nEY = nE;
    const int nos = 4;
    LM.reserve(nEX*nEY);
    int e = 0;
    int chi = nE - 1;

    for(int col = 0; col < nEX; col++) {
        for(int row = 0; row < nEY; row++) {
            e = col + row*nEX;

            std::vector<int> global_p;
            global_p.reserve(nos);

            for(int k = 0; k < nos; k++) {
                
                if ( k <= 1) {
                    global_p.push_back(k + e + row + 1);
                } else {
                    global_p.push_back(k + e + row + chi + 1);
                }
            }
            LM.push_back(global_p);
        }
    }

    return LM;
}

matrix build_Ce2D(std::function<double(double,double)> f) {
    const int size = 4;
    double* Ce = new double[size*size];

    for(int row = 0; row < size; row++ ) {
            for(int col = 0; col < size; col++) {
                Ce[row + col*size] = gauss_legendre2D([&row,&col, f](double x,double y){return f(x,y)*n1pol2D(row,x,y)*n1pol2D(col,x,y);},5);
            }
    }

    matrix result = {Ce, size, size};
    return result;
}


matrix build_Ae2D(std::function<double(double,double)> f) {
    const int size = 4;
    double* Ae = new double[size*size];

    for(int row = 0; row < size; row++ ) {
            for(int col = 0; col < size; col++) {
                Ae[row + col*size] = gauss_legendre2D([&row,&col, f](double x,double y){return f(x,y)*n1dx_pol2D(row,x,y)*n1dx_pol2D(col,x,y)+f(x,y)*n1dy_pol2D(row,x,y)*n1dy_pol2D(col,x,y);},5);
            }
    }

    matrix result = {Ae, size, size};
    return result;
}



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
                Be[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n3pol(row,x)*n3dpol(col,x);},4);
                break;

            case 2: // Quadrático
                Be[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n2pol(row,x)*n2dpol(col,x);},2);
                break;
            
            default: // Linear
                Be[row + col*width] = gauss_legendre([&row,&col, f](double x){return f(x)*n1pol(row,x)*n1dpol(col,x);},2);
            }
        }
    }

    matrix result = {Be, size, size};


    return result;
}

vector build_Fe(int n, double a, double b, std::function<double(double)> f) {

    assert(n >= 1);

    int size = n + 1;

    double* Fe = new double[size];

    for(int row = 0; row < size; row++) {
        switch (n)
            {
            case 3: // Cúbico
                Fe[row] = gauss_legendre([&row, &a, &b, f](double e){return f( (b-a)/2*(e+1)+a )*n3pol(row,e);},3);
                break;

            case 2: // Quadrático
                Fe[row] = gauss_legendre([&row, &a, &b, f](double e){return f( (b-a)/2*(e+1)+a )*n2pol(row,e);},3);
                break;
            
            default: // Linear
                Fe[row] = gauss_legendre([&row, &a, &b, f](double e){return f( (b-a)/2*(e+1)+a )*n1pol(row,e);},3);
            }
    }

    vector result = {Fe, size};

    return result;
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


matrix build_K(std::vector<double> X,
               std::vector<std::vector<int>> LM, 
               std::function<double(double)> alpha, 
               std::function<double(double)> beta, 
               std::function<double(double)> gamma, 
               int n) {
    int N = X.size();

    // Número de pontos para K
    // Ne = N-1 (número de elementos)
    // np = n+1 (número de pontos por elementos)
    // Número de pontos = Ne*np -(Ne-1)
    // pois tem que retirar os pontos que conectam os elementos
    int size_K = (N-1)*(n+1) - (N-2);

    double he = X[1]-X[0];

    auto Ae = build_Ae(n,alpha);
    auto Be = build_Be(n,beta);
    auto Ce = build_Ce(n,gamma);

    double* K = new double[size_K*size_K];
    for(int i = 0; i < size_K*size_K; i++) K[i] = 0;

    int e_index = 0 ; // index do elemento
    int local_i, local_j;
    for(auto e: LM) {

        for(int local_i = 0; local_i < e.size(); local_i ++) {
            for(int local_j = 0; local_j < e.size(); local_j ++) {
                auto global_i = e[local_i];
                auto global_j = e[local_j];

                K[global_i + global_j*size_K] += 2/he*Ae.m[local_i + local_j*Ae.nrows] + he/2*Ce.m[local_i + local_j*Ce.nrows] + Be.m[local_i + local_j*Be.nrows];
            }
        }
    }

    K[0] += 2/he*Ae.m[Ae.nrows*Ae.ncols-1] + he/2*Ce.m[Ce.nrows*Ce.ncols-1] + Be.m[Be.ncols*Be.nrows-1];
    K[size_K*size_K-1] += 2/he*Ae.m[0] + he/2*Ce.m[0] + Be.m[0];


    matrix result = {K, size_K, size_K};

    delete[] Ae.m;
    delete[] Be.m;
    delete[] Ce.m;

    return result;
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

    
    

    double* F = new double[sizeF];
    
    for(int i = 0; i < sizeF; i++) F[i] = 0;

    int e_index = 0 ; // index do elemento
    int local_i;
    for(auto e: LM) {

        auto Fe = build_Fe(n, X[e.front()], X[e.back()], f);

        for(auto global_i: e) {
            local_i = global_i - n*e_index;

            F[global_i] += he/2*Fe.m[local_i];
        }
        e_index++;

        delete[] Fe.m;
    }

    auto Fe = build_Fe(n, X[LM.front().front()], X[LM.front().back()], f);
    F[0] += he/2*Fe.m[Fe.n-1];
    delete[] Fe.m;

    Fe = build_Fe(n, X[LM.back().front()], X[LM.back().back()], f);
    F[sizeF-1] += he/2*Fe.m[0];
    delete[] Fe.m;

    vector result = {F, sizeF};

    

    return result;
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
