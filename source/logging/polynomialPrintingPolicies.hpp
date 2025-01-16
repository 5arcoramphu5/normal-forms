#pragma once

#include "../containers/Polynomial.hpp"

template<typename T>
concept PolynomialPrintingPolicy = requires (Polynomial<double> p){
    { T::polyToString(p) } -> std::same_as<std::string>;
};

const std::string defaultVars[] = {"x1", "x2", "x3", "x4"};
template<ArithmeticType Coeff>
std::string toString(Polynomial<Coeff> polynomial, const std::string vars[] = defaultVars);

struct SymbolicPolynomialPrinting {

    template<ArithmeticType Coeff>
    static std::string polyToString(const Polynomial<Coeff> &p)
    { return toString(p); }
};

template<ArithmeticType Coeff>
std::string toCoefficientString(Polynomial<Coeff> polynomial);

struct CoefficientPolynomialPrinting {

    template<ArithmeticType Coeff>
    static std::string polyToString(const Polynomial<Coeff> &p)
    { return toCoefficientString(p); }
};

#include "polynomialPrintingPolicies.tpp"