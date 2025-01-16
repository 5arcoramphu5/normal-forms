#include "capd/capdlib.h"
#include "helperFunctions.hpp"
#include "../logging/logging.hpp"

using namespace capd;
using namespace std;

Polynomial<Complex> getTaylorSeries(const CMap &function, int degree)
{
    CJet taylor(function.imageDimension(), function.dimension(), degree);
    function(CVector({0, 0, 0, 0}), taylor);
    return taylor;
}

Polynomial<Complex> projP(const Polynomial<Complex> &poly, int upToDegree)
{
    int maxDeg = upToDegree != -1 ? std::min((int)poly.degree(), upToDegree) : poly.degree();
    Polynomial<Complex> result(4, 4, maxDeg);

    for(int i = 0; i < maxDeg; ++i)
        for(int j = 0; 2*i+2*j < maxDeg ; ++j)
        {   
            Multiindex i1({i+1, i, j, j});
            result(i1)[0] = poly(0, i1); // P_1

            Multiindex i2({i, i+1, j, j});
            result(i2)[1] = poly(1, i2); // P_2

            Multiindex i3({i, i, j+1, j});
            result(i3)[2] = poly(2, i3); // P_3

            Multiindex i4({i, i, j, j+1});
            result(i4)[3] = poly(3, i4); // P_4
        }

    return result;
}

Polynomial<Complex> projR(const Polynomial<Complex> &poly, int upToDegree)
{
    int maxDeg = upToDegree != -1 ? std::min((int)poly.degree(), upToDegree) : poly.degree();
    Polynomial<Complex> result(4, 4, maxDeg);

    for(int deg = 0; deg <= maxDeg; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            for(int i = 0; i < 4; ++i)
            {
                if(
                    (i == 0 && index[0]-1 == index[1] && index[2] == index[3]) ||
                    (i == 1 && index[0] == index[1]-1 && index[2] == index[3]) ||
                    (i == 2 && index[0] == index[1] && index[2]-1 == index[3]) ||
                    (i == 3 && index[0] == index[1] && index[2] == index[3]-1) )
                    continue;
                    
                result(index)[i] = poly(i, index);
            }
        }while(index.hasNext());
    }

    return result;
}

CVector gamma(int p, int q, Complex lambda1, Complex lambda2)
{
    return CVector({
        Complex(p-1, 0)*lambda1 + Complex(q, 0)*lambda2,
        Complex(p+1, 0)*lambda1 + Complex(q, 0)*lambda2,
        Complex(p, 0)*lambda1 + Complex(q-1, 0)*lambda2,
        Complex(p, 0)*lambda1 + Complex(q+1, 0)*lambda2
    });
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

PairMap<Polynomial<Complex>> pqCoefficients(const Polynomial<Complex> & poly, int upToDegree)
{
    PairMap<Polynomial<Complex>> coefficients;

    for(int deg = 0; deg <= upToDegree; ++deg)
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

                if (!coefficients.contains(pair_pq))
                {
                    Polynomial<Complex> newElem(4, 2, upToDegree);
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

Polynomial<Complex> operatorL(const Polynomial<Complex> Psi, const Polynomial<Complex> &N, const CMatrix &lambda)
{
    return D(Psi) * N - (Polynomial<Complex>)(lambda*Psi);
}

PolynomialMatrix<Complex, 4> D(const Polynomial<Complex> &F)
{
    PolynomialMatrix<Complex, 4> result(F.degree());

    const Multiindex zero(F.dimension());

    for(int i = 0; i < F.imageDimension(); ++i)
        for(int j = 0; j < F.dimension(); ++j)
        {   
            // d f_i / d x_j
            for(int deg = 1; deg <= F.degree(); ++deg)
            {
                Multiindex index(zero);
                index[0] = deg;
                do
                {
                    if(index[j] > 0)
                    {
                        Multiindex newIndex(index);
                        newIndex[j] -= 1;
                        result[i][j](0, newIndex) = Complex(index[j], 0) * F(i, index);
                    }

                }while(index.hasNext());
            }
        }
    
    return result;
}