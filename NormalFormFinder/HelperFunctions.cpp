#include "capd/capdlib.h"
#include "HelperFunctions.h"

using namespace capd;
using namespace std;

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, PolynomialOf4Variables4 *reminder)
{    
    if(taylor.dimension() != 4)
        throw new runtime_error("Taylor series dimension is invalid.");

    Complex linearArr[4][4];
    Multiindex index({1, 0, 0, 0});
    int i = 0;
    do
    {
        for(int j = 0; j < 4; ++j)
            linearArr[i][j] = taylor(index)[j];
        i++;
    }while(index.hasNext());

    *linearPart = CMatrix(linearArr);

    for(int deg = 0; deg <= taylor.degree(); ++deg)
    {
        if(deg != 1)
        {
            Multiindex index({deg, 0, 0, 0});
            do
            {
                for(int j = 0; j < 4; ++j)
                    (*reminder)(j).setCoeff(index, taylor(index)[j]);
            }while(index.hasNext());
        }
    }
}

CJet getTaylorSeries(const CMap &function, int degree)
{
    CJet taylor(function.imageDimension(), function.dimension(), degree);
    function(CVector({0, 0, 0, 0}), taylor);
    return taylor;
}

CVector getEigenvalues(const DMatrix &matrix)
{
    auto dim = matrix.dimension();
    if(dim.first != dim.second)
        throw runtime_error("Matrix is not square.");

    DVector eigenValuesR(dim.first), eigenValuesI(dim.first);
    CVector eigenValues(dim.first);
    alglib::computeEigenvalues(matrix, eigenValuesR, eigenValuesI);
    for(int i = 0; i < 4; ++i)
    {
        eigenValues[i] = Complex(eigenValuesR[i], eigenValuesI[i]);
    }

    return eigenValues;
}

PolynomialOf4Variables4 proj_P(const PolynomialOf4Variables4 &poly)
{
    int deg = poly.getDegree();
    PolynomialOf4Variables4 result(deg);

    for(int i = 0; i < deg; ++i)
        for(int j = 0; 2*i+2*j < deg ; ++j)
        {   
            if(j > 0)
            {
                result(0).setCoeff(i+1, i, j, j, poly(0).getCoeff(i+1, i, j, j)); // P_1
                result(1).setCoeff(i, i+1, j, j, poly(1).getCoeff(i, i+1, j, j)); // P_2
            }

            if(i > 0)
            {
                result(2).setCoeff(i, i, j+1, j, poly(2).getCoeff(i, i, j+1, j)); // P_3
                result(3).setCoeff(i, i, j, j+1, poly(3).getCoeff(i, i, j, j+1)); // P_4
            }
        }

    return result;
}

PolynomialOf4Variables4 proj_R(const PolynomialOf4Variables4 &poly)
{
    PolynomialOf4Variables4 result(poly.getDegree());
    auto proj = proj_P(poly);
    for(int i = 0; i < 4; ++i)
        result(i) = poly(i) - proj(i);
    return result;
}
