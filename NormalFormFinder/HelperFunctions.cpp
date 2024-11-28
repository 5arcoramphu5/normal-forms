#include "capd/capdlib.h"
#include "HelperFunctions.h"

using namespace capd;
using namespace std;

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, CJet *reminder)
{    
    if(taylor.dimension() != 4)
        throw new runtime_error("Taylor series dimension is invalid.");

    *reminder = taylor;

    Complex linearArr[4][4];
    for(int i = 0; i < 4; ++i)
    {
        int indexArr[] = {0, 0, 0, 0};
        indexArr[i] = 1;
        Multiindex index(4, indexArr);

        for(int j = 0; j < 4; ++j)
            linearArr[i][j] = taylor(j, index);

        (*reminder)(index) = {0, 0, 0, 0};
    }

    *linearPart = CMatrix(linearArr);
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

CJet projP(const CJet &poly)
{
    int deg = poly.degree();
    CJet result(4, deg);

    for(int i = 0; i < deg; ++i)
        for(int j = 0; 2*i+2*j < deg ; ++j)
        {   
            if(j > 0)
            {
                Multiindex i1({i+1, i, j, j});
                result(i1)[0] = poly(0, i1); // P_1

                Multiindex i2({i, i+1, j, j});
                result(i2)[1] = poly(1, i2); // P_2
            }

            if(i > 0)
            {
                Multiindex i3({i, i, j+1, j});
                result(i3)[2] = poly(2, i3); // P_3

                Multiindex i4({i, i, j, j+1});
                result(i4)[3] = poly(3, i4); // P_4
            }
        }

    return result;
}

CJet projR(const CJet &poly)
{
    CJet result(4, poly.degree());
    auto proj = projP(poly);
    return poly - proj;
}

string toString(CJet polynomial, string var1, string var2, string var3, string var4)
{    
    const Multiindex zero({0, 0, 0, 0});
    string result = "";

    for(int i = 0; i < polynomial.dimension(); ++i)
    {
        stringstream ss;
        for(int deg = 0; deg <= polynomial.degree(); ++deg)
        {
            Multiindex index({deg, 0, 0, 0});
            do
            {
                Complex coeff = polynomial(i, index);
                if(coeff != Complex(0, 0))
                {
                    if(index != zero)
                    {
                        if(coeff != Complex(1, 0))
                        {
                            if(coeff == Complex(-1, 0)) ss << "-";
                            else ss << coeff;
                        }
                    }
                    else ss << coeff;

                    if(index[0] > 0) ss << var1<<(index[0] == 1 ? "" : "^"+to_string(index[0]));
                    if(index[1] > 0) ss << var2<<(index[1] == 1 ? "" : "^"+to_string(index[1]));
                    if(index[2] > 0) ss << var3<<(index[2] == 1 ? "" : "^"+to_string(index[2]));
                    if(index[3] > 0) ss << var4<<(index[3] == 1 ? "" : "^"+to_string(index[3]));

                    ss << " + ";
                }
            }while(index.hasNext());
        }

        string str = ss.str();
        if(str.length() == 0)
            str = "0";
        else 
            str = str.substr(0, str.length()-3);

        result += str + "\n";
    }

    return result;
}
