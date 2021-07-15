#pragma once


inline double n1pol_a(double x) {
    return 0.5-x/2.0;
}

inline double n1pol_b(double x) {
    return 0.5+x/2.0;
}

double n1pol(int i, double x);
double n1dpol(int i, double x);

double n2pol(int i, double x);
double n2dpol(int i, double x);

double n3pol(int i, double x);
double n3dpol(int i, double x);

inline double map_x_to_e(double x, double a, double b);

double pol(int order, int i, double x);
double dpol(int order, int i, double x);

void interpolate(double x_eval, double a, double b, int order, double* W, double &u, double &du);


//
//
//     2______________3
//     |              |
//     |              |
//     |              |
//     |              |
//     |______________|
//     0              1 
//
//

double inline pol2d_a(double x, double y) {
    return n1pol_a(x)*n1pol_a(y);
}

double inline pol2d_b(double x, double y) {
    return n1pol_b(x)*n1pol_a(y);
}

double inline pol2d_c(double x, double y) {
    return n1pol_a(x)*n1pol_b(y);
}

double inline pol2d_d(double x, double y) {
    return n1pol_b(x)*n1pol_b(y);
}

double n1pol2D(int i, double x, double y);


inline double n1dpol_a(double x) {
    return -0.5;
}

inline double n1dpol_b(double x) {
    return 0.5;
}

double inline dx_pol2d_a(double x, double y) {
    return n1dpol_a(x)*n1pol_a(y);
}

double inline dy_pol2d_a(double x, double y) {
    return n1pol_a(x)*n1dpol_a(y);
}

double inline dx_pol2d_b(double x, double y) {
    return n1dpol_b(x)*n1pol_a(y);
}

double inline dy_pol2d_b(double x, double y) {
    return n1pol_b(x)*n1dpol_a(y);
}

double inline dx_pol2d_c(double x, double y) {
    return n1dpol_a(x)*n1pol_b(y);
}

double inline dy_pol2d_c(double x, double y) {
    return n1pol_a(x)*n1dpol_b(y);
}

double inline dx_pol2d_d(double x, double y) {
    return n1dpol_b(x)*n1pol_b(y);
}

double inline dy_pol2d_d(double x, double y) {
    return n1pol_b(x)*n1dpol_b(y);
}

double n1dx_pol2D(int i, double x, double y);
double n1dy_pol2D(int i, double x, double y);