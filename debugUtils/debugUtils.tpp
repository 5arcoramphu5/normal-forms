#pragma once

#include "debugUtils.hpp"

// print in format compatible with Mathematica
// ostream& operator<<(ostream& os, const Complex& c)
// {
//     os << "(" << c.real() << " + " << c.imag() << "I)";
//     return os;
// }

template<ArithmeticType Coeff>
std::string toString(Polynomial<Coeff> polynomial, const std::string vars[])
{    
    const capd::Multiindex zero(polynomial.dimension());
    std::string result = "\t";

    for(int i = 0; i < polynomial.imageDimension(); ++i)
    {
        std::stringstream ss;
        for(int deg = 0; deg <= polynomial.degree(); ++deg)
        {
            capd::Multiindex index(polynomial.dimension());
            index[0] = deg;

            do
            {
                Coeff coeff = polynomial(i, index);
                if(coeff != (Coeff)0)
                {
                    if(index != zero)
                    {
                        if(coeff != (Coeff)1)
                        {
                            if(coeff == (Coeff)-1) ss << "-";
                            else ss << coeff << " ";
                        }
                    }
                    else ss << coeff << " ";
                    
                    for(int j = 0; j < 4; ++j)
                        if(index[j] > 0) ss << vars[j] << (index[j] == 1 ? "" : "^"+std::to_string(index[j])) << " ";

                    ss << " + ";
                }
            }while(index.hasNext());
        }

        std::string str = ss.str();
        if(str.length() == 0)
            str = "0";
        else 
            str = str.substr(0, str.length()-3);

        result += str + "\n\t";
    }

    return result;
}

template<ArithmeticType Coeff>
std::string toCoefficientString(Polynomial<Coeff> polynomial)
{
    std::stringstream result;
    result << "{\n";
    for(int i = 0; i < polynomial.imageDimension(); ++i)
    {
        for(int deg = 0; deg <= polynomial.degree(); ++deg)
        {
            capd::Multiindex index(polynomial.dimension());
            index[0] = deg;

            do
            {
                Coeff coeff = polynomial(i, index);
                result << coeff << " ";
            }while(index.hasNext());

            result << "\n";
        }
        if(i != polynomial.imageDimension()-1)
            result << ",\n";
    }
    result << "}";

    return result.str();
}