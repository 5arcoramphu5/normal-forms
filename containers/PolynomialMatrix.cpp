#include "PolynomialMatrix.h"

#include "capd/capdlib.h"

using namespace std;
using namespace capd;

void multiplyAndAdd(Polynomial &accumulator, int accumulatorIndex, const Polynomial &p1, int p1Index, const Polynomial &p2, int p2Index)
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

Polynomial operator*(const PolynomialMatrix<4> &jetMatrix, const Polynomial &jet)
{
    if(jet.dimension() != 4) throw runtime_error("invalid dimension for polynomial matrix multiplication.");

    Polynomial result(4, 4, jetMatrix.degree() + jet.degree());

    for(int row = 0; row < 4; ++row)
        for(int i = 0; i < 4; ++i)
            multiplyAndAdd(result, row, jetMatrix[row][i], 0, jet, i);

    return result;
}