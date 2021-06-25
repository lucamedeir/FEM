#include<iostream>
#include<vector>
//#include <gtest/gtest.h>
#include <args-parser/all.hpp>
#include "builders.hpp"
#include "utils.hpp"
#include "interpol.hpp"
#include <cmath>

inline double dsolution(double x) {
    return 1.0/2.0-x;
}

inline double solution(double x) {
    return -x*(x-1.0)/2.0;
}

double errorL2(vector u, std::vector<std::vector<int>> LM) {

    auto newX = linspace(0.0, 1.0, u.n);
    double I2 = 0;
    double r,dr;
    double rs,drs;
    double a,b;


    int order = 0;
    double eval,deval;
    for(auto e: LM) {
        order = e.size()-1;

        a = newX[e[0]];
        b = newX[e[order]];

        auto W = u.m+e[0];

        eval = 0;
        deval = 0;
        for(auto gauss: g_roots[3]) {

            r = 0;
            dr = 0;

            for(int i = 0; i <= order; i++) {
                r += W[i] * pol(order,i,gauss.abscissa);
                dr += W[i] * dpol(order,i,gauss.abscissa);
            }

            dr *= 2/(b-a);

            rs = solution((b-a)/2*(gauss.abscissa+1)+a);
            drs = dsolution((b-a)/2*(gauss.abscissa+1)+a);
            
            eval += gauss.weight*( (rs - r)*(rs - r) );
            deval += gauss.weight*( (drs - dr)*(drs - dr) );
        }

        I2 += eval*(b-a)/2 + deval*(b-a)/2;
    }

    return sqrt(I2);
}

void execute_and_print_error(const int NumberOfPoints,
                            const int order,
                            const double alpha,
                            const double beta,
                            const double gamma,
                            const double u_initial, 
                            const double u_final) {

    auto X = linspace(0.0, 1.0, NumberOfPoints);
    auto LM = build_LM(NumberOfPoints, order);

    auto alpha_f = [&alpha](double x) {return alpha;};
    auto beta_f = [&beta](double){return beta;};
    auto gamma_f = [&gamma](double x) {return gamma;};
    auto f = [](double x) {return x;};

    auto K = build_K(X, LM, alpha_f, beta_f, gamma_f, order);
    auto F = build_F(X, LM, f, order);

    set_boundary(K,F,order,order,u_initial,u_final);

    auto u = solve(K,F);

    std::cout << errorL2(u, LM) << std::endl;

    delete[] K.m;
    delete[] F.m;
    delete[] u.m;

}

void execute(const int NumberOfPoints,
            const int order,
            const double alpha,
            const double beta,
            const double gamma,
            const double u_initial, 
            const double u_final) {

    auto X = linspace(0.0, 1.0, NumberOfPoints);
    auto LM = build_LM(NumberOfPoints, order);


    auto alpha_f = [&alpha](double x) {return alpha;};
    auto beta_f = [&beta](double){return beta;};
    auto gamma_f = [&gamma](double x) {return gamma;};
    auto f = [](double x) {return x;};

    auto K = build_K(X, LM, alpha_f, beta_f, gamma_f, order);
    auto F = build_F(X, LM, f, order);

    set_boundary(K,F,order,order,u_initial,u_final);

    auto u = solve(K,F);
    print(u,LM);


    

    delete[] K.m;
    delete[] F.m;
    delete[] u.m;

}

void execute_and_interpolate(const int NumberOfPoints,
                            const int order,
                            const double alpha,
                            const double beta,
                            const double gamma,
                            const double u_initial, 
                            const double u_final,
                            const int Nint) {

    auto X = linspace(0.0, 1.0, NumberOfPoints);
    auto LM = build_LM(NumberOfPoints, order);

    auto alpha_f = [&alpha](double x) {return alpha;};
    auto beta_f = [&beta](double){return beta;};
    auto gamma_f = [&gamma](double x) {return gamma;};
    auto f = [](double x) {return x;};

    auto K = build_K(X, LM, alpha_f, beta_f, gamma_f, order);
    auto F = build_F(X, LM, f, order);
    

    set_boundary(K,F,order,order,u_initial,u_final);

    auto u = solve(K,F);

    auto newX = linspace(0.0, 1.0, u.n);

    interpolate_and_print(newX, u, LM, Nint);

    delete[] K.m;
    delete[] F.m;
    delete[] u.m;

}


int parser(int argc, char** argv, int &N, int &order, double &alpha, double &beta, double &gamma, double &ui, double &uf, int &Nint, bool &inter, bool &printError) {
    try {
        Args::CmdLine cmd(argc, argv);

        Args::Arg NumberOfPointsArg( SL('n'), SL("npoints"), true, true );

        NumberOfPointsArg.setDescription( SL("Number of Points. Any natural number bigger or equal to 2."));
        NumberOfPointsArg.setLongDescription( SL("Number of Points. Number of nodes in the domain. Any natural number bigger than or equal to 2. "
          "The Number of elements will be equal to this minus one. When value is 2, therefore there is only one element."));
    
        Args::Arg AlphaArg( SL('a'), SL("alpha"), true, true);
        AlphaArg.setDescription( SL("Constant value that multiplies Ae matrix."));
        AlphaArg.setLongDescription( SL("Constant value that multiplies Ae matrix, the u'v' dx part." ));

        Args::Arg BetaArg( SL('b'), SL("beta"), true, true);
        AlphaArg.setDescription( SL("Constant value that multiplies Be matrix."));
        AlphaArg.setLongDescription( SL("Constant value that multiplies Be matrix, the u v' dx part." ));

        Args::Arg GammaArg( SL('c'), SL("gamma"), true, true);
        GammaArg.setDescription( SL("Constant value that multiplies Ce matrix."));
        GammaArg.setLongDescription( SL("Constant value that multiplies Ce matrix, the u v dx part." ));

        Args::Arg UInitialArg( SL('i'), SL("initial"), true, false);
        UInitialArg.setDescription( SL("Initial value of solution. Default is 0."));
        UInitialArg.setLongDescription( SL("Initial value of solution. Default is 0. Dirichlet Boundary problem." ));

        Args::Arg UFinalArg( SL('f'), SL("final"), true, false);
        UFinalArg.setDescription( SL("Final value of solution. Default is 0."));
        UFinalArg.setLongDescription( SL("Final value of solution. Default is 0. Dirichlet Boundary problem." ));

        Args::Arg OrderArg( SL('o'),SL("order"), true, false);
        OrderArg.setDescription( SL("Default order for all elements. Default is one."));
        OrderArg.setLongDescription (SL("Interpolation order for all elements. Default is one. Order of 1 means linear, 2 means quadratic and 3 means cubic."));

        Args::Arg InterpolateArg( SL('t'), SL("interpolate"), true, false);
        InterpolateArg.setDescription( SL("Number of points for the interpolation. If set, the program returns the x, u and du."));
        InterpolateArg.setLongDescription( SL("Number of points for the interpolation. If set, the program returns the abscissa vector, the evaluation of the solution at the given abscissa points as well as the derivate."));

        Args::Arg ErrorL2Arg( SL('e'), SL("errorl2"), false, false);
        ErrorL2Arg.setDescription( SL("Return the L2 error."));
        ErrorL2Arg.setLongDescription( SL("Return the L2 error."));

        Args::Help help;
        help.setExecutable(argv[0]);
        help.setAppDescription(SL("This application solves the poisson problem with finite element methods on the domain [0,1]. Returns a table with weights of every element and its interpolation order."));

        cmd.addArg(NumberOfPointsArg);
        cmd.addArg(AlphaArg);
        cmd.addArg(BetaArg);
        cmd.addArg(GammaArg);
        cmd.addArg(OrderArg);
        cmd.addArg(UInitialArg);
        cmd.addArg(UFinalArg);
        cmd.addArg(InterpolateArg);
        cmd.addArg(ErrorL2Arg);
        cmd.addArg(help);

        cmd.parse();

        N = std::stoi(NumberOfPointsArg.value());
        alpha = std::stod(AlphaArg.value());
        beta = std::stod(BetaArg.value());
        gamma = std::stod(GammaArg.value());

        if( OrderArg.isDefined() )
            order = std::stoi(OrderArg.value());
        if( UInitialArg.isDefined() )
            ui = std::stod(UInitialArg.value());
        if( UFinalArg.isDefined() )
            uf = std::stod(UFinalArg.value());

        inter = InterpolateArg.isDefined();
        printError = ErrorL2Arg.isDefined();

        if (inter) 
            Nint = std::stoi(InterpolateArg.value());

    }
    catch( const Args::HelpHasBeenPrintedException & )
	{
        return 2;
	}
	catch( const Args::BaseException & x )
	{
		Args::outStream() << x.desc() << SL( "\n" );
        return 1;
	}

    return 0;
}


int main(int argc, char **argv) {/*
    ::testing::InitGoogleTest(&argc, argv);
    auto test = RUN_ALL_TESTS();*/


    int NumberOfPoints = 2;
    int NInt = 4;
    bool printInterpolation = false;
    bool printError = true;
    int order = 1;
    double alpha = 1;
    double beta = 0;
    double gamma = 0;
    double ui = 0;
    double uf = 0;
   
    int result = parser(argc, argv, NumberOfPoints, order, alpha, beta, gamma, ui, uf, NInt, printInterpolation, printError);
    //int result = 0;


    if(result == 0) {
        if (printInterpolation) {
            execute_and_interpolate(NumberOfPoints, order, alpha, beta, gamma, ui, uf, NInt);
        } else if (printError) {
            execute_and_print_error(NumberOfPoints, order, alpha, beta, gamma, ui, uf);
        } else {
            execute(NumberOfPoints, order, alpha, beta, gamma, ui, uf);
        }
    }


    //return test;
    return 0;
}