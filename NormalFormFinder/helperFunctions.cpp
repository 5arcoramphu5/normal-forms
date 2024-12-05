#include "capd/capdlib.h"
#include "helperFunctions.h"

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

CJet projP(const CJet &poly, int upToDeg)
{
    int maxDeg = upToDeg != -1 ? std::min((int)poly.degree(), upToDeg) : poly.degree();
    CJet result(4, maxDeg);

    for(int i = 0; i < maxDeg; ++i)
        for(int j = 0; 2*i+2*j < maxDeg ; ++j)
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

CJet projR(const CJet &poly, int upToDeg)
{
    int maxDeg = upToDeg != -1 ? std::min((int)poly.degree(), upToDeg) : poly.degree();
    CJet result(4, maxDeg);

    for(int deg = 0; deg <= maxDeg; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            for(int i = 0; i < 4; ++i)
            {
                bool skip = false;
                switch(i)
                {
                    case 0:
                        if(index[0]+1 == index[1] && index[2] == index[3])
                            skip = true;
                        break;
                    case 1:
                        if(index[0] == index[1]+1 && index[2] == index[3])
                            skip = true;
                        break;
                    case 2:
                        if(index[0] == index[1] && index[2]+1 == index[3])
                            skip = true;
                        break;
                    case 3:
                        if(index[0] == index[1] && index[2] == index[3]+1)
                            skip = true;
                        break;
                }

                if(!skip)
                    result(index)[i] = poly(i, index);
            }
        }while(index.hasNext());
    }

    return result;
}

vector<Complex> gamma(int p, int q, Complex lambda1, Complex lambda2)
{
    vector<Complex> data(4);

    data[0] = Complex(p-1, 0)*lambda1 + Complex(q, 0)*lambda2;
    data[1] = Complex(p+1, 0)*lambda1 + Complex(q, 0)*lambda2;
    data[2] = Complex(p, 0)*lambda1 + Complex(q-1, 0)*lambda2;
    data[3] = Complex(p, 0)*lambda1 + Complex(q+1, 0)*lambda2;
    
    return data;
}

bool isNonzero(CColumnVector columnVector)
{
    for(Complex x : columnVector)
    {
        if(x.real() != 0 || x.imag() != 0) 
            return true;
    }
    return false;
}

std::unordered_map<std::pair<int,int>,CJet,hash_pair> pqCoefficients(const CJet & poly)
{
    unordered_map<pair<int, int>, CJet, hash_pair> coefficients;
    int degree = poly.degree();

    for(int deg = 0; deg <= degree; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            auto coeffVector = poly(index);

            if(isNonzero(coeffVector))
            {
                int p = index[0] - index[1]; // j - k
                int q = index[2] - index[3]; // l - m
                auto pair_pq = make_pair(p, q);

                if (coefficients.find(pair_pq) == coefficients.end()) // not present in dictionary
                {
                    CJet newElem(4, 2, 2*degree);
                    coefficients.insert(make_pair(pair_pq, newElem));
                }

                Multiindex coefficientIndex({index[1], index[3]}); // (k, m)
                coefficients[pair_pq](coefficientIndex) += coeffVector;
            }
        }
        while(index.hasNext());
    }

    return coefficients;
}

