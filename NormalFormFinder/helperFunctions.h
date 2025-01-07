#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include "../containers/Polynomial.h"
#include "../containers/PolynomialMatrix.h"

#include <unordered_map>

Polynomial getTaylorSeries(const CMap &function, int degree);

Polynomial projP(const Polynomial &poly, int upToDegree = -1);

Polynomial projR(const Polynomial &poly, int upToDegree = -1);

CVector gamma(int p, int q, capd::Complex lambda1, capd::Complex lambda2);

// required to create unordered_map with keys of type pair<int, int>
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        size_t hash1 = std::hash<T1>{}(p.first);
        size_t hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
};

std::unordered_map<std::pair<int, int>, Polynomial, hash_pair> pqCoefficients(const Polynomial &poly, int upToDegree);

Polynomial operatorL(const Polynomial Psi, const Polynomial &N, const CMatrix &lambda);

PolynomialMatrix<4> D(const Polynomial &F);

#endif