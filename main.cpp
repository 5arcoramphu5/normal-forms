#include <iostream>
#include "capd/capdlib.h"
#include "NormalFormFinder/NormalFormFinder.h" 
#include "NormalFormFinder/PseudoNormalForm.h"
#include "Polynomials/PolynomialOf4Variables.h"

using namespace std;
using namespace capd;

int main()
{
    PolynomialOf4Variables poly(2, 5, 0);
    poly.extend(1);
    poly.set_coeff(1, 1, 1, 0, 3);
    cout << poly.toString("x", "y", "z", "u") <<endl;


    return 0;
}