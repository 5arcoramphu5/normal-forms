#pragma once

#include "polynomialPrintingPolicies.hpp"
#include "capd/fields/lib.h"

template<ArithmeticType T>
T round(T value, double precisionMargin)
{
    if(std::abs(value) <= precisionMargin)
        return 0;
    return value;
}

template<ArithmeticType Coeff>
std::string toString(Polynomial<Coeff> polynomial, const std::string vars[], double precisionMargin)
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
                coeff = round(coeff, precisionMargin);
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
std::string toCoefficientString(Polynomial<Coeff> polynomial, double precisionMargin)
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
                coeff = round(coeff, precisionMargin);
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