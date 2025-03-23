#pragma once 

#include "Polynomial.hpp"

template<ArithmeticType Coeff>
Polynomial<Coeff> Polynomial<Coeff>::fromToDegree(int degreeFrom, int degreeTo) const
{
    if(degreeFrom < 0 || degreeTo < 0 || degreeFrom > degreeTo) throw std::runtime_error("invalid degrees");

    if(degreeTo > this->degree())
        degreeTo = this->degree();
    
    Polynomial<Coeff> result(this->imageDimension(), this->dimension(), degreeTo);

    for(int deg = degreeFrom; deg <= degreeTo; ++deg)
    {
        capd::Multiindex index(this->dimension());
        index[0] = deg;

        do
        {
            for(int i = 0; i < this->imageDimension(); ++i)
                result(i, index) = (*this)(i, index);
        }while(index.hasNext());
    }

    return result;
}

template<ArithmeticType Coeff>
Polynomial<Coeff> operator+(const Polynomial<Coeff> &p1, const Polynomial<Coeff> &p2)
{
    int resultDegree = std::max(p1.degree(), p2.degree());
    Polynomial<Coeff> result(4, 4, resultDegree);

    for(int deg = 0; deg <= resultDegree; ++deg)
    {
        capd::Multiindex index({deg, 0, 0, 0});
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

template<ArithmeticType Coeff>
Polynomial<Coeff> operator-(const Polynomial<Coeff> &p1, const Polynomial<Coeff> &p2)
{
    int resultDegree = std::max(p1.degree(), p2.degree());
    Polynomial<Coeff> result(4, 4, resultDegree);

    for(int deg = 0; deg <= resultDegree; ++deg)
    {
        capd::Multiindex index({deg, 0, 0, 0});
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

template<ArithmeticType Coeff>
Coeff compositionProduct(const Polynomial<Coeff> &second, const capd::Multiindex& mi, const capd::Multipointer& a, int p, int k)
{
    Coeff result = 0;
    const auto is = capd::Multipointer::generateList(p,k);

    auto e = is.end();
    for(auto b = is.begin(); b != e; ++b)
    {
        auto bt = b->begin(), et = b->end();
        auto ib = mi.begin();

        capd::Multipointer delta = a.subMultipointer(*bt);

        if(delta.dimension() > second.degree())
            continue;

        Coeff temp = second(*ib,delta) * Coeff(delta.factorial(), 0);
        ++bt;
        ++ib;

        for( ; bt != et; ++bt)
        {
            capd::Multipointer delta = a.subMultipointer(*bt);

            if(delta.dimension() > second.degree())
            {
                temp = 0;
                break;
            }

            temp *= second(*ib,delta) * Coeff(delta.factorial(), 0);
            ++ib;
        }

        result += temp;
    }

    return result;
}

template<ArithmeticType Coeff>
void composition(const Polynomial<Coeff> &first, const Polynomial<Coeff> &second, Polynomial<Coeff> &result, const capd::Multipointer& a)
{
    typename capd::Multiindex::IndicesSet listIndices;
    capd::Multiindex::generateList(result.dimension(), result.degree(), listIndices);

    int p = a.module();

    result(a).clear();
    int minK = 1;
    for(int k = minK; k <= p; ++k)
    {
        auto e = listIndices[k-1].end();
        for(auto b = listIndices[k-1].begin(); b != e; ++b)
        {
            capd::Multipointer mp(b->dimension(),b->begin());

            std::sort(mp.begin(),mp.end());
            auto product = compositionProduct(second, *b, a, p, k);
            result(a) += first(mp) * product * (Coeff)mp.factorial();
        }
    }

    result(a) /= (Coeff)a.factorial();
}

// modified version of substitutionPowerSeries(...), should work for second.degree() < first.degree()
template<ArithmeticType Coeff>
void polynomialComposition(const Polynomial<Coeff> &first, const Polynomial<Coeff> &second, Polynomial<Coeff> &result)
{
    for(unsigned i=1;i<=first.degree();++i)
    {
        capd::Multipointer a = first.first(i);
        do
        {
            composition(first, second, result, a);
        }
        while(first.hasNext(a));
    }
    result() = first();
}

template<ArithmeticType Coeff>
Polynomial<Coeff> taylorSeriesAtPoint(const Map<Coeff> &map, const Vector<Coeff> &point, int degree)
{
    if(degree == 0) // segmentation fault in CAPD in this case
    {
        Polynomial<Coeff> result(map.imageDimension(), map.dimension(), 0);
        result(capd::Multiindex(map.imageDimension())) = map(point);
        return result;
    }

    Polynomial<Coeff> result(map.imageDimension(), map.dimension(), degree);
    map(point, result);
    return result;
}

// Taylor series expansion of division of two polynomials C^2 -> C^4, filling only coefficients of degree pased as argument
template<ArithmeticType Coeff>
Polynomial<Coeff> polynomialDivision(const Polynomial<Coeff> &numerator, const Polynomial<Coeff> &denominator, int degree)
{
    // Taylor expansion of x/y function
    Map<Coeff> division(
        "var: x, y; "
        "fun: x/y;",
        degree);

    int polyDeg = numerator.degree();
    Polynomial<Coeff> result(4, 2, degree);

    for(int i = 0; i < 4; ++i)
    {
        Polynomial<Coeff> p(2, 2, polyDeg);
        for(int deg = 0; deg <= polyDeg; ++deg)
        {
            capd::Multiindex index({deg, 0});
            do
            {
                p(0, index) = numerator(i, index);
                p(1, index) = denominator(i, index);
            }while(index.hasNext());
        }

        capd::Multiindex zero({0, 0});
        auto divisionTaylor = taylorSeriesAtPoint(division, {numerator(i, zero), denominator(i, zero)}, degree);

        Polynomial<Coeff> composition(1, 2, degree);
        polynomialComposition(divisionTaylor, p, composition);

        capd::Multiindex index({degree, 0});
        do
        {
            result(i, index) = composition(0, index);
        }while(index.hasNext());
    }

    return result;
}

template<ArithmeticType Coeff>
Polynomial<Coeff> toPolynomial(const capd::vectalg::Matrix<Coeff, 0, 0> &linearPart, const capd::vectalg::Vector<Coeff, 0> &constant, int degree)
{
    auto dim = linearPart.dimension();
    if(constant.dimension() != dim.first)
        throw std::runtime_error("invalid dimensions of matrix and vector");

    Polynomial<Coeff> result(dim.first, dim.second, 1);
    
    for(int i = 0; i < dim.first; ++i)
    {
        for(int j = 0; j < dim.second; ++j)
        {
            capd::Multiindex index(dim.second);
            index[j] = 1;

            result(i, index) = linearPart[i][j];
        }

        result(i, capd::Multiindex(dim.second)) = constant[i];
    }

    return result;
}