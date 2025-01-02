#include <iostream>

#include "typedefs.h"
#include "debugUtils/debugUtils.h"
#include "NormalFormFinder/NormalFormFinder.cpp"
using namespace std;
using namespace capd;

#define MAX_DERIVATIVE 10
#define METHOD_DEGREE 2

void diagonal_matrix_test()
{
    CMap f(
        "par: p1, p2, p3, p4, a1, a2, a3;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "p1*x1 + a1*x1*x3*x4,"
            "p2*x2 + a2*x2*x4,"
            "p3*x3 + a3*x2*x2,"
            "p4*x4 + x1*x1*x1 + x2*x3;", 
        MAX_DERIVATIVE);

    f.setParameter("p1", Complex(1, 1));
    f.setParameter("p2", Complex(-1, -1));
    f.setParameter("p3", Complex(1, -1));
    f.setParameter("p4", Complex(-1, 1));
    f.setParameter("a1", 2);
    f.setParameter("a2", 5);
    f.setParameter("a3", 1i);

    CVector p({1, 2, 3, 4});
    CMatrix lambda({ {Complex(1, 1), 0, 0, 0}, {0, Complex(-1, -1), 0, 0}, {0, 0, Complex(1, -1), 0}, {0, 0, 0, Complex(-1, 1)} });
    CMatrix J({ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} });
    CMatrix invJ({ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} });

    NormalFormFinder<Logger<VerbosityLevel::Diagnostic>> finder(METHOD_DEGREE, f, p, lambda, J, invJ);
    PseudoNormalForm normalForm = finder.calculatePseudoNormalForm();
}

void henon_heiles_test()
{
    CMap f(
        "par: a1, lambda;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "x2,"
            "-x1 - a1*lambda*x1*x3,"
            "x4,"
            "-x3 - lambda*(x1^2 - x3^2);", 
        MAX_DERIVATIVE);

    f.setParameter("a1", 2);
    f.setParameter("lambda", 1);

    CMatrix lambda({ {Complex(0, -1), 0, 0, 0}, {0, Complex(0, -1), 0, 0}, {0, 0, Complex(0, 1), 0}, {0, 0, 0, Complex(0, 1)} });
    CMatrix J({ {0, 0, Complex(0, -0.5), Complex(0.5, 0)}, {Complex(0, -0.5), Complex(0.5, 0), 0, 0}, {0, 0, Complex(0, 0.5), Complex(0.5, 0)}, {Complex(0, 0.5), Complex(0.5, 0), 0, 0} });
    CMatrix invJ({ {0, Complex(0, 1), 0, Complex(0, -1)}, {0, 1, 0, 1}, {Complex(0, 1), 0, Complex(0, -1), 0}, {0, 1, 0, 1} });

    CVector p({0, 0, 0, 0});

    NormalFormFinder<Logger<VerbosityLevel::Diagnostic>> finder(METHOD_DEGREE, f, p, lambda, J, invJ);
    PseudoNormalForm normalForm = finder.calculatePseudoNormalForm();
}

int main()
{
    diagonal_matrix_test();
    henon_heiles_test();

    return 0;
}