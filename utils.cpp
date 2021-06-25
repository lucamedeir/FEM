
#include "utils.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <algorithm>
#include <fstream>
#include "interpol.hpp"



std::vector<std::vector<gaussian_point>> g_roots = {
    { // n = 2
        {1.0,-0.5773502691896257},
        {1.0,0.5773502691896257}
    },
    { // n = 3
        {0.5555555555555556,-0.7745966692414834},
        {0.8888888888888888,0.0},
        {0.5555555555555556,0.7745966692414834}
    },
    { // n = 4
        {0.3478548451374538, -0.8611363115940526},
        {0.6521451548625461, -0.3399810435848563},
        {0.6521451548625461, 0.3399810435848563},
        {0.3478548451374538, 0.8611363115940526}
    },
    { // n = 5
        {0.2369268850561891, -0.9061798459386640},
        {0.4786286704993665, -0.5384693101056831},
        {0.5688888888888889, 0.0000000000000000},
        {0.4786286704993665, 0.5384693101056831},
        {0.2369268850561891, 0.9061798459386640}
    }
};



void print( vector u, std::vector<std::vector<int>> LM) {
    std::cout << "element,order,w1,w2,w3,w4," << std::endl;

    int e = 0;
    for(auto global_p: LM) {
        std::cout << e << ",";
        
        std::cout << global_p.size()-1 << ",";
        for(auto p : global_p) {
            std::cout << u.m[p] << ",";
        }
        std::cout << std::endl;
        e++;
    }

}

void save( vector u, std::vector<std::vector<int>> LM, const char* filename) {
    std::ofstream savefile;
    savefile.open (filename);
    savefile << "element,order,w1,w2,w3,w4," << std::endl;

    int e = 0;
    for(auto global_p: LM) {
        savefile << e << ",";
        
        savefile << global_p.size()-1 << ",";
        for(auto p : global_p) {
            savefile << u.m[p] << ",";
        }
        savefile << std::endl;
        e++;
    }

    savefile.close();
}

double gauss_legendre(std::function<double(double)> f, int n) {
    double eval = 0;

    assert(n>=2);

    for(auto gauss: g_roots[n-2]) {
        eval += gauss.weight*f(gauss.abscissa);
    } 

    return eval;
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

std::vector<double> linspace(double initial_value, double final_value, int N) {
    std::vector<double> X(N);
    double Dx = 1.0/(N-1);
    std::generate(X.begin(), X.end(), [n = 0, &Dx]() mutable { return n++ * Dx; } );
    return X;
}

TEST(BuildingTest, linspaceTest) {

    auto X = linspace(0.0, 1.0, 5);

    std::vector<double> XTest = {0.0, 0.25, 0.5, 0.75, 1.0};

    int i = 0;
    for(auto v: X) {
        ASSERT_EQ(v,XTest[i++]) << "position i in vector X:" << i;
    }
}

void evaluate_solution(double x, std::vector<double> X, vector u, std::vector<std::vector<int>> LM, double &r, double &dr) {
    int order = 1;
    int element_index = 0;
    double a =0;
    double b = 0;

    for(auto e: LM) {
        order = e.size() - 1;

        if( x>=  X[element_index*order] && x <= X[order*(element_index +1)] ) {
            a = X[element_index*order];
            b = X[order*(element_index +1)];
            interpolate(x, a, b, order, u.m + e[0], r, dr);
            break;
        } 
        element_index++;
    }
}

void interpolate_and_print(std::vector<double> X, vector u, std::vector<std::vector<int>> LM, int Nint) {
    auto Y = linspace(0,1,Nint);
    double r;
    double dr;

    std::cout << "x,u,du" << std::endl;
    for(auto y: Y) {
        evaluate_solution(y,X,u,LM,r,dr);
        std::cout << y << "," << r << "," << dr << std::endl;
    }
}