#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "../typedefs.h"

// a wrapper for CJet, enabling operations on Jets of different degrees
class Polynomial : public CJet
{
    public:
        using CJet::Jet;

        Polynomial(const CJet &other) : CJet(other) {}

        Polynomial fromToDegree(int degreeFrom, int degreeTo) const;

        Polynomial inline reminderPart() const
        { return fromToDegree(2, degree()); }
};

Polynomial operator+(const Polynomial &p1, const Polynomial &p2);
Polynomial operator-(const Polynomial &p1, const Polynomial &p2);

void polynomialComposition(const Polynomial &first, const Polynomial &second, Polynomial& result);

Polynomial polynomialDivision(const Polynomial &numerator, const Polynomial &denominator, int degree);

Polynomial toPolynomial(const CMatrix &linearPart, const CVector &constant, int degree);

#endif