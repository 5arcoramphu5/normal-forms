#include "PolynomialMatrix.h"

#include "capd/capdlib.h"

using namespace std;
using namespace capd;

template<ArithmeticType Coeff>
void multiplyAndAdd(Polynomial<Coeff> &accumulator, int accumulatorIndex, const Polynomial<Coeff> &p1, int p1Index, const Polynomial<Coeff> &p2, int p2Index)
{
    for(int deg1 = 0; deg1 <= p1.degree(); ++deg1)
    {
        Multiindex index1({deg1, 0, 0, 0});
        do
        {
            for(int deg2 = 0; deg2 <= p2.degree(); ++deg2)
            {
                Multiindex index2({deg2, 0, 0, 0});
                do
                {
                    Multiindex indexSum({index1[0]+index2[0], index1[1]+index2[1], index1[2]+index2[2], index1[3]+index2[3]});
                    accumulator(accumulatorIndex, indexSum) += p1(p1Index, index1) * p2(p2Index, index2); 
                }while(index2.hasNext());
            }
        }while(index1.hasNext());
    }
}

template<ArithmeticType Coeff, int N>
Polynomial<Coeff> operator*(const PolynomialMatrix<Coeff, N> &jetMatrix, const Polynomial<Coeff> &jet)
{
    if(jet.dimension() != N) throw runtime_error("invalid dimension for polynomial matrix multiplication.");

    Polynomial<Coeff> result(N, N, jetMatrix.degree() + jet.degree());

    for(int row = 0; row < N; ++row)
        for(int i = 0; i < N; ++i)
            multiplyAndAdd(result, row, jetMatrix[row][i], 0, jet, i);

    return result;
}