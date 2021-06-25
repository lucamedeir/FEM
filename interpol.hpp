#pragma once

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