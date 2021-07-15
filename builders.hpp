#pragma once
#include <functional>

typedef struct matrix {
    double* m;
    int nrows;
    int ncols;
} matrix;

typedef struct vector {
    double* m;
    int n;
} vector;

matrix build_Ce2D(std::function<double(double,double)> f);
matrix build_Ae2D(std::function<double(double,double)> f);


matrix build_Ce(int n, std::function<double(double)> f);
matrix build_Ae(int n, std::function<double(double)> f);
matrix build_Be(int n, std::function<double(double)> f);

vector build_Fe(int n, double a, double b, std::function<double(double)> f);
std::vector<std::vector<int>> build_LM(int nX, int n);
matrix build_K(std::vector<double> X, std::vector<std::vector<int>> LM, std::function<double(double)> alpha, std::function<double(double)> beta, std::function<double(double)> gamma, int n);
vector build_F(std::vector<double> X, std::vector<std::vector<int>> LM, std::function<double(double)> f,  int n);
vector solve(matrix K, vector F);
void set_boundary(matrix K,vector F, double ni, double nf,  double ui, double uf);