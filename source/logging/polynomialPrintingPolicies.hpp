#pragma once

#include "../containers/Polynomial.hpp"

template<typename T>
concept PolynomialPrintingPolicy = requires (Polynomial<double> p, double precisionMargin){
    { T::polyToString(p, precisionMargin) } -> std::same_as<std::string>;
};

const std::string defaultVars[] = {"x1", "x2", "x3", "x4"};
template<ArithmeticType Coeff>
std::string toString(Polynomial<Coeff> polynomial, const std::string vars[] = defaultVars, double precisionMargin = 0);

struct SymbolicPolynomialPrinting {

    template<ArithmeticType Coeff>
    static std::string polyToString(const Polynomial<Coeff> &p, double precisionMargin = 0)
    { return toString(p, defaultVars, precisionMargin);}
};

struct CoefficientPolynomialPrinting {

    template<ArithmeticType Coeff>
    static std::string polyToString(const Polynomial<Coeff> &p, double precisionMargin = 0);
};

// only for capd::fields::Complex<T>
struct MathematicaFormatPolynomialPrinting {

    template<ArithmeticType Coeff>
    static std::string polyToString(const Polynomial<Coeff> &p, double precisionMargin = 0);
};

#include "polynomialPrintingPolicies.tpp"