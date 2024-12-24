#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include <unordered_map>
#include <vector>

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, CJet *reminder);

CJet getTaylorSeries(const CMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

CJet projP(const CJet &poly, int upToDegree = -1);

CJet projR(const CJet &poly, int upToDegree = -1);

std::vector<capd::Complex> gamma(int p, int q, capd::Complex lambda1, capd::Complex lambda2);

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

std::unordered_map<std::pair<int, int>, CJet, hash_pair> pqCoefficients(const CJet &poly, int upToDegree);

CJet polyDivision(const CJet &numerator, const CJet &denominator);

CJet fromToDegree(const CJet &poly, int degreeFrom, int degreeTo);

CJet operatorL(const CJet Psi, const CJet &N, const CMatrix &lambda);

CJet jetSubstraction(const CJet &p1, const CJet &p2);

CJet jetAddition(const CJet &p1, const CJet &p2);

template<int N>
class CJetMatrix : public std::array<std::array<CJet, N>, N>
{
    int matrixDegree;

    public:
    CJetMatrix(int degree) : matrixDegree(degree)
    { 
        for(auto &array : *this)
            for(auto &jet : array)
                jet = CJet(1, N, degree);
    }

    int degree() const { return matrixDegree; }
};

CJetMatrix<4> D(const CJet &F);

CJet inline reminderPart(const CJet &poly)
{ return fromToDegree(poly, 2, poly.degree()); }

#endif