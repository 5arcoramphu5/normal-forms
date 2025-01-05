#include "capd/capdlib.h"
#include "helperFunctions.h"

using namespace capd;
using namespace std;

CJet getTaylorSeries(const CMap &function, int degree)
{
    CJet taylor(function.imageDimension(), function.dimension(), degree);
    function(CVector({0, 0, 0, 0}), taylor);
    return taylor;
}

CJet projP(const CJet &poly, int upToDegree)
{
    int maxDeg = upToDegree != -1 ? std::min((int)poly.degree(), upToDegree) : poly.degree();
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

CJet projR(const CJet &poly, int upToDegree)
{
    int maxDeg = upToDegree != -1 ? std::min((int)poly.degree(), upToDegree) : poly.degree();
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
                        if(index[0]-1 == index[1] && index[2] == index[3])
                            skip = true;
                        break;
                    case 1:
                        if(index[0] == index[1]-1 && index[2] == index[3])
                            skip = true;
                        break;
                    case 2:
                        if(index[0] == index[1] && index[2]-1 == index[3])
                            skip = true;
                        break;
                    case 3:
                        if(index[0] == index[1] && index[2] == index[3]-1)
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

std::unordered_map<std::pair<int,int>,CJet,hash_pair> pqCoefficients(const CJet & poly, int upToDegree)
{
    unordered_map<pair<int, int>, CJet, hash_pair> coefficients;

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

                if (coefficients.find(pair_pq) == coefficients.end()) // not present in dictionary
                {
                    CJet newElem(4, 2, upToDegree);
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

// division of two polynomials C^2 -> C^4, only passed degree coefficients
CJet polyDivision(const CJet &numerator, const CJet &denominator, int degree)
{
    // TODO: implement proper division
    CJet result(numerator);
    Multiindex zero({0, 0});

    Multiindex index({degree, 0});
        do
        {
            for(int i = 0; i < 4; ++i)
            {
                if(result(i, index) == Complex(0, 0)) continue;
                result(i, index) /= denominator(i, zero);
            }

        }while(index.hasNext());
    
    return result;
}

CJet fromToDegree(const CJet &poly, int degreeFrom, int degreeTo)
{
    if(degreeFrom < 0 || degreeTo < 0 || degreeFrom > degreeTo) throw runtime_error("invalid degrees");

    if(degreeTo > poly.degree())
        degreeTo = poly.degree();
    
    CJet result(poly.imageDimension(), poly.dimension(), degreeTo);

    for(int deg = degreeFrom; deg <= degreeTo; ++deg)
    {
        Multiindex index(poly.dimension());
        index[0] = deg;

        do
        {
            for(int i = 0; i < poly.imageDimension(); ++i)
                result(i, index) = poly(i, index);
        }while(index.hasNext());
    }

    return result;
}

void multiplyAndAdd(CJet &accumulator, int accumulatorIndex, const CJet &p1, int p1Index, const CJet &p2, int p2Index)
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

CJet multiply(const CJetMatrix<4> &jetMatrix, const CJet &jet)
{
    if(jet.dimension() != 4) throw runtime_error("invalid dimension for jet matrix multiplication.");

    CJet result(4, 4, jetMatrix.degree() + jet.degree());

    for(int row = 0; row < 4; ++row)
        for(int i = 0; i < 4; ++i)
            multiplyAndAdd(result, row, jetMatrix[row][i], 0, jet, i);

    return result;
}

CJet operatorL(const CJet Psi, const CJet &N, const CMatrix &lambda)
{
    return jetSubstraction(multiply(D(Psi), N), lambda*Psi);
}

CJet jetSubstraction(const CJet &p1, const CJet &p2) // TODO: to be deleted
{
    int resultDegree = std::max(p1.degree(), p2.degree());
    CJet result(4, 4, resultDegree);

    for(int deg = 0; deg <= resultDegree; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            for(int i = 0; i < 4; ++i)
            {
                if(deg <= p1.degree())
                    result(i, index) += p1(i, index);
                
                if(deg <= p2.degree())
                    result(i, index) -= p2(i, index);
            }
        }while(index.hasNext());
    }
    return result;
}

CJet jetAddition(const CJet &p1, const CJet &p2) // TODO: to be deleted
{
    int resultDegree = std::max(p1.degree(), p2.degree());
    CJet result(4, 4, resultDegree);

    for(int deg = 0; deg <= resultDegree; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            for(int i = 0; i < 4; ++i)
            {
                if(deg <= p1.degree())
                    result(i, index) += p1(i, index);
                
                if(deg <= p2.degree())
                    result(i, index) += p2(i, index);
            }
        }while(index.hasNext());
    }
    return result;
}

CJet jetAddition(const CMatrix &linearPart, const CVector &constant, int degree) // to be deleted
{
    auto dim = linearPart.dimension();
    if(constant.dimension() != dim.first)
        throw runtime_error("invalid dimensions of matrix and vector");

    CJet result(dim.first, dim.second, /* debug */ degree);
    
    for(int i = 0; i < dim.first; ++i)
    {
        for(int j = 0; j < dim.second; ++j)
        {
            Multiindex index(dim.second);
            index[j] = 1;

            result(i, index) = linearPart[i][j];
        }

        result(i, Multiindex(dim.second)) = constant[i];
    }

    return result;
}

CJetMatrix<4> D(const CJet &F)
{
    CJetMatrix<4> result(F.degree());

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
