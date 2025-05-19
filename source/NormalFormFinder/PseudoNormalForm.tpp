#pragma once
#include "PseudoNormalForm.hpp"

#include <iostream>

CVector PseudoNormalForm::solution(double t, CVector initialPoint) const
{
    CVector constPart({initialPoint[0]*initialPoint[1], initialPoint[2]*initialPoint[3]}); // value of xi*eta, mu*nu - constant with respect to t
    auto a1_const = a1(constPart)[0];
    auto a2_const = a2(constPart)[0];
    return CVector(
        {
            exp(a1_const*t) * initialPoint[0],
            exp(-a1_const*t) * initialPoint[1],
            exp(a2_const*t) * initialPoint[2],
            exp(-a2_const*t) * initialPoint[3]
        }
    );
}

void PseudoNormalForm::serialize(std::ostream &stream) const
{
    phi.serialize(stream);
    n.serialize(stream);
    b.serialize(stream);
    a1.serialize(stream);
    a2.serialize(stream);
}

PseudoNormalForm PseudoNormalForm::deserialize(std::istream &stream)
{
    auto phi = Polynomial<capd::Complex>::deserialize(stream);
    auto n = Polynomial<capd::Complex>::deserialize(stream);
    auto b = Polynomial<capd::Complex>::deserialize(stream);
    auto a1 = Polynomial<capd::Complex>::deserialize(stream);
    auto a2 = Polynomial<capd::Complex>::deserialize(stream);
    return PseudoNormalForm(phi, n, b, a1, a2);
}