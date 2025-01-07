#include "Polynomial.h"

#include "capd/capdlib.h"

using namespace std;
using namespace capd;

Polynomial Polynomial::fromToDegree(int degreeFrom, int degreeTo) const
{
    if(degreeFrom < 0 || degreeTo < 0 || degreeFrom > degreeTo) throw runtime_error("invalid degrees");

    if(degreeTo > degree())
        degreeTo = degree();
    
    Polynomial result(imageDimension(), dimension(), degreeTo);

    for(int deg = degreeFrom; deg <= degreeTo; ++deg)
    {
        Multiindex index(dimension());
        index[0] = deg;

        do
        {
            for(int i = 0; i < imageDimension(); ++i)
                result(i, index) = (*this)(i, index);
        }while(index.hasNext());
    }

    return result;
}

Polynomial operator+(const Polynomial &p1, const Polynomial &p2)
{
    int resultDegree = std::max(p1.degree(), p2.degree());
    Polynomial result(4, 4, resultDegree);

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

Polynomial operator-(const Polynomial &p1, const Polynomial &p2)
{
    int resultDegree = std::max(p1.degree(), p2.degree());
    Polynomial result(4, 4, resultDegree);

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

Complex compositionProduct(const Polynomial &second, const Multiindex& mi, const Multipointer& a, int p, int k)
{
    Complex result = 0;
    const auto is = Multipointer::generateList(p,k);

    auto e = is.end();
    for(auto b = is.begin(); b != e; ++b)
    {
        auto bt = b->begin(), et = b->end();
        auto ib = mi.begin();

        Multipointer delta = a.subMultipointer(*bt);

        if(delta.dimension() > second.degree())
            continue;

        Complex temp = second(*ib,delta) * Complex(delta.factorial(), 0);
        ++bt;
        ++ib;

        for( ; bt != et; ++bt)
        {
            Multipointer delta = a.subMultipointer(*bt);

            if(delta.dimension() > second.degree())
            {
                temp = 0;
                break;
            }

            temp *= second(*ib,delta) * Complex(delta.factorial(), 0);
            ++ib;
        }

        result += temp;
    }

    return result;
}

void composition(const Polynomial &first, const Polynomial &second, Polynomial &result, const Multipointer& a)
{
    typename Multiindex::IndicesSet listIndices;
    Multiindex::generateList(result.dimension(), result.degree(), listIndices);

    int p = a.module();

    result(a).clear();
    int minK = 1;
    for(int k = minK; k <= p; ++k)
    {
        auto e = listIndices[k-1].end();
        for(auto b = listIndices[k-1].begin(); b != e; ++b)
        {
            Multipointer mp(b->dimension(),b->begin());

            sort(mp.begin(),mp.end());
            auto product = compositionProduct(second, *b, a, p, k);
            result(a) += first(mp) * product *  Complex(mp.factorial(), 0);
        }
    }

    result(a) /= Complex(a.factorial(), 0);
}

// modified version of substitutionPowerSeries(...), should work for second.degree() < first.degree()
void polynomialComposition(const Polynomial &first, const Polynomial &second, Polynomial& result)
{
    for(unsigned i=1;i<=first.degree();++i)
    {
        Multipointer a = first.first(i);
        do
        {
            composition(first, second, result, a);
        }
        while(first.hasNext(a));
    }
    result() = first();
}

// division of two polynomials C^2 -> C^4, only passed degree coefficients
Polynomial polynomialDivision(const Polynomial &numerator, const Polynomial &denominator, int degree)
{
    // TODO: implement proper division
    Polynomial result(numerator);
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

Polynomial toPolynomial(const CMatrix &linearPart, const CVector &constant, int degree) // to be deleted
{
    auto dim = linearPart.dimension();
    if(constant.dimension() != dim.first)
        throw runtime_error("invalid dimensions of matrix and vector");

    Polynomial result(dim.first, dim.second, 1);
    
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