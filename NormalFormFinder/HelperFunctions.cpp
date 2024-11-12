#include "capd/capdlib.h"
#include "HelperFunctions.h"

using namespace capd;
using namespace std;

void getLinearPartWithReminder(const DJet &taylor, DMatrix *linearPart, PolynomialOf4Variables4 *reminder)
{    
    if(taylor.dimension() != 4)
        throw new runtime_error("Taylor series dimension is invalid.");

    double linearArr[4][4];
    Multiindex index({1, 0, 0, 0});
    int i = 0;
    do
    {
        for(int j = 0; j < 4; ++j)
            linearArr[i][j] = taylor(index)[j];
        i++;
    }while(index.hasNext());

    *linearPart = DMatrix(linearArr);

    for(int deg = 0; deg <= taylor.degree(); ++deg)
    {
        if(deg != 1)
        {
            Multiindex index({deg, 0, 0, 0});
            do
            {
                for(int j = 0; j < 4; ++j)
                    (*reminder)(j).set_coeff(index, taylor(index)[j]);
            }while(index.hasNext());
        }
    }
}

DJet getTaylorSeries(const DMap &function, int degree)
{
    DJet taylor(function.imageDimension(), function.dimension(), degree);
    function(DVector({0, 0}), taylor);
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