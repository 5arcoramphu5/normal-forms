#include <iostream>

#include "typedefs.h"
#include "debugUtils/debugUtils.h"
#include "NormalFormFinder/NormalFormFinder.cpp"
using namespace std;
using namespace capd;

int main()
{
    CVector x({1.5, 2, 1, 0.75});
    int maxDerivativeOrder = 5;
    CMap f(
        "par: p1, p2, p3, p4, a1, a2, a3;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "p1*x1 + a1*x1*x3*x4,"
            "p2*x2 + a2*x2*x4,"
            "p3*x3 + a3*x2*x2,"
            "p4*x4 + x1*x1*x1 + x2*x3;", 
        maxDerivativeOrder);

    f.setParameter("p1", Complex(1, 1));
    f.setParameter("p2", Complex(-1, -1));
    f.setParameter("p3", Complex(1, -1));
    f.setParameter("p4", Complex(-1, 1));
    f.setParameter("a1", 2);
    f.setParameter("a2", 5);
    f.setParameter("a3", 1i);

    NormalFormFinder<Logger<VerbosityLevel::Diagnostic>> finder(2, f, x);
    PseudoNormalForm normalForm = finder.calculatePseudoNormalForm();

    // cout << "Phi:\n" << toString(normalForm.getPhi()) << endl;
    // cout << "N:\n" << toString(normalForm.getN()) << endl;
    // cout << "B:\n" << toString(normalForm.getB()) << endl;
    // checkPseudoNormalCondition(normalForm);

    return 0;
}