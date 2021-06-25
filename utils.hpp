#pragma once

#include<functional>
#include<vector>
#include"builders.hpp"


typedef struct gaussian_point {
    double weight;
    double abscissa;
} gaussian_point;

extern std::vector<std::vector<gaussian_point>> g_roots;



double gauss_legendre(std::function<double(double)> f, int n);
std::vector<double> linspace(double initial_value, double final_value, int N);
void save( vector u, std::vector<std::vector<int>> LM, const char* filename);
void print( vector u, std::vector<std::vector<int>> LM);
void interpolate_and_print(std::vector<double> X, vector u, std::vector<std::vector<int>> LM, int Nint);
void evaluate_solution(double x, std::vector<double> X, vector u, std::vector<std::vector<int>> LM, double &r, double &dr);