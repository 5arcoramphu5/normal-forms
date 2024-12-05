#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include <unordered_map>
#include <vector>

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, CJet *reminder);

CJet getTaylorSeries(const CMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

CJet projP(const CJet &poly, int upToDeg = -1);

CJet projR(const CJet &poly, int upToDeg = -1);

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

std::unordered_map<std::pair<int, int>, CJet, hash_pair> pqCoefficients(const CJet &poly);

#endif