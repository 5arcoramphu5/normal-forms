#pragma once

#include <unordered_map>

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

template<typename T>
class PairMap : public std::unordered_map<std::pair<int, int>, T, hash_pair>
{
    public:
        inline bool contains(const std::pair<int, int> &pair)
        { return this->find(pair) != this->end(); }
};