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

template <ArithmeticType Coeff>
void Polynomial<Coeff>::serialize(std::ostream &stream) const
{
    stream << this->degree() << " " << this->dimension() << " " << this->imageDimension() << "\n";
    for(int i = 0; i < this->imageDimension(); ++i)
    {
        for(int deg = 0; deg <= this->degree(); ++deg)
        {
            capd::Multiindex index(this->dimension());
            index[0] = deg;
            do{
                stream << (*this)(i, index) << " ";
            }
            while(index.hasNext());
            stream << "\n";
        }
    }
}

template <ArithmeticType Coeff>
Polynomial<Coeff> Polynomial<Coeff>::deserialize(std::istream &stream)
{
    int degree, dim, imDim;
    stream >> degree >> dim >> imDim;
    Polynomial<Coeff> poly(imDim, dim, degree);

    for(int i = 0; i < imDim; ++i)
    {
        for(int deg = 0; deg <= degree; ++deg)
        {
            capd::Multiindex index(dim);
            index[0] = deg;
            do{
                stream >> poly(i, index);
            }
            while(index.hasNext());
        }
    }
    return poly;
}


template<ArithmeticType Coeff>
Polynomial<Coeff> operator+(const Polynomial<Coeff> &p1, const Polynomial<Coeff> &p2)
{
    if(p1.dimension() != p2.dimension() || p1.imageDimension() != p2.imageDimension())
        throw std::runtime_error("invalid polynomial dimensions");

    int resultDegree = std::max(p1.degree(), p2.degree());
    Polynomial<Coeff> result(4, 4, resultDegree);

    for(int deg = 0; deg <= resultDegree; ++deg)
    {
        capd::Multiindex index(p1.dimension());
        index[0] = deg;
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
    if(p1.dimension() != p2.dimension() || p1.imageDimension() != p2.imageDimension())
        throw std::runtime_error("invalid polynomial dimensions");

    int resultDegree = std::max(p1.degree(), p2.degree());
    Polynomial<Coeff> result(p1.imageDimension(), p1.dimension(), resultDegree);

    for(int deg = 0; deg <= resultDegree; ++deg)
    {
        capd::Multiindex index(p1.dimension());
        index[0] = deg;
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

// Taylor series expansion of division of two polynomials C^2 -> C^4, filling only coefficients of degree pased as argument
template<ArithmeticType Coeff>
Polynomial<Coeff> polynomialDivision(const Polynomial<Coeff> &numerator, const Polynomial<Coeff> &denominator, int degree)
{
    int polyDeg = numerator.degree();

    // Taylor expansion of x/y function
    Map<Coeff> division(
        "var: x, y; "
        "fun: x/y;",
        polyDeg);

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

        Polynomial<Coeff> composition = division(p);

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