#include "interpol.hpp"


inline double n1pol_a(double x) {
    return 0.5-x/2.0;
}

inline double n1pol_b(double x) {
    return 0.5+x/2.0;
}

double n1pol(int i, double x) {
    switch(i) {
        case 0: return n1pol_a(x);
    }

    return n1pol_b(x);
}

inline double n1dpol_a(double x) {
    return -0.5;
}

inline double n1dpol_b(double x) {
    return 0.5;
}

double n1dpol(int i, double x) {
    switch(i) {
        case 0: return n1dpol_a(x);
    }

    return n1dpol_b(x);
}

inline double n2pol_a(double x) {
    return x*x/2.0 -x/2.0;
}

inline double n2dpol_a(double x) {
    return x-1.0/2.0;
}

inline double n2pol_b(double x) {
    return 1.0-x*x;
}

inline double n2dpol_b(double x) {
    return -2.0*x;
}

inline double n2pol_c(double x) {
    return x*x/2.0 + x/2.0;
}

inline double n2dpol_c(double x) {
    return x+1.0/2.0;
}

double n2pol(int i, double x) {
    switch(i) {
        case 0: return n2pol_a(x);
        case 1: return n2pol_b(x);
    }

    return n2pol_c(x);
}

double n2dpol(int i, double x) {
    switch(i) {
        case 0: return n2dpol_a(x);
        case 1: return n2dpol_b(x);
    }

    return n2dpol_c(x);
}

inline double n3pol_a(double x) {
    return -0.5625*x*x*x + 0.5625*x*x + 0.0625*x - 0.0625;
}

inline double n3dpol_a(double x) {
    return -3.0*0.5625*x*x + 2.0*0.5625*x + 0.0625;
}

inline double n3pol_b(double x) {
    return 1.6875*x*x*x - 0.5625*x*x - 1.6875*x + 0.5625;
}

inline double n3dpol_b(double x) {
    return 3.0*1.6875*x*x - 2.0*0.5625*x - 1.6875;
}

inline double n3pol_c(double x) {
    return -1.6875*x*x*x - 0.5625*x*x + 1.6875*x + 0.5625;
}

inline double n3dpol_c(double x) {
    return -3.0*1.6875*x*x - 2.0*0.5625*x + 1.6875;
}

inline double n3pol_d(double x) {
    return 0.5625*x*x*x + 0.5625*x*x - 0.0625*x - 0.0625;
}

inline double n3dpol_d(double x) {
    return 3.0*0.5625*x*x + 2.0*0.5625*x - 0.0625;
}

double n3pol(int i, double x) {
    switch(i) {
        case 0: return n3pol_a(x);
        case 1: return n3pol_b(x);
        case 2: return n3pol_c(x);
    }

    return n3pol_d(x);
}

double n3dpol(int i, double x) {
    switch(i) {
        case 0: return n3dpol_a(x);
        case 1: return n3dpol_b(x);
        case 2: return n3dpol_c(x);
    }

    return n3dpol_d(x);
}

double pol(int order, int i, double x) {
    switch(order) {
        case 3: return n3pol(i, x);
        case 2: return n2pol(i, x);
    }

    return n1pol(i,x);
}

double dpol(int order, int i, double x) {
    switch(order) {
        case 3: return n3dpol(i, x);
        case 2: return n2dpol(i, x);
    }

    return n1dpol(i,x);
}

inline double map_x_to_e(double x, double a, double b) {
    return 2*(x-a)/(b-a)-1;
}

void interpolate(double x_eval, double a, double b, int order, double* W, double &u, double &du) {
    u = 0;
    du = 0;
    for(int i = 0; i <= order; i++) {
        u += W[i] * pol(order,i,map_x_to_e(x_eval,a,b));
        du += W[i] * dpol(order,i,map_x_to_e(x_eval,a,b));
    }
    du *= 2/(b-a);
}