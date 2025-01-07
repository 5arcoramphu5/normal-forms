#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "capd/vectalg/lib.h"
#include "../typedefs.h"

template<typename T>
concept ArithmeticType = requires (T a, T b) {
    a + b; a += b;
    a - b; a -= b;
    a * b; a *= b;
    a / b; a /= b;
};

// a wrapper for CJet, enabling operations on Jets of different degrees
template<ArithmeticType Coeff>
class Polynomial : public capd::diffAlgebra::Jet<
                            capd::vectalg::Matrix<Coeff, 0, 0>, 
                            0>
{
    public:
        using MatrixType = capd::vectalg::Matrix<Coeff, 0, 0>;
        using JetType = capd::diffAlgebra::Jet<MatrixType, 0>;

        using JetType::JetType;

        Polynomial(const JetType &other) : JetType(other) {}

        Polynomial<Coeff> fromToDegree(int degreeFrom, int degreeTo) const;

        Polynomial<Coeff> inline reminderPart() const
        { return fromToDegree(2, this->degree()); }
};

template<ArithmeticType Coeff>
Polynomial<Coeff> operator+(const Polynomial<Coeff> &p1, const Polynomial<Coeff> &p2);

template<ArithmeticType Coeff>
Polynomial<Coeff> operator-(const Polynomial<Coeff> &p1, const Polynomial<Coeff> &p2);

template<ArithmeticType Coeff>
void polynomialComposition(const Polynomial<Coeff> &first, const Polynomial<Coeff> &second, Polynomial<Coeff> &result);

template<ArithmeticType Coeff>
Polynomial<Coeff> polynomialDivision(const Polynomial<Coeff> &numerator, const Polynomial<Coeff> &denominator, int degree);

template<ArithmeticType Coeff>
Polynomial<Coeff> toPolynomial(const capd::vectalg::Matrix<Coeff, 0, 0> &linearPart, const capd::vectalg::Vector<Coeff, 0> &constant, int degree);

#endif