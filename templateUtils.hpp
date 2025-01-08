#pragma once

#include <string>

template <typename T>
concept Streamable = requires(std::ostream &os, T value) {
    { os << value } -> std::convertible_to<std::ostream&>;
};

template <typename T>
concept Invokable = requires(T t) {
    t();
};

template<typename T>
concept ArithmeticType = requires (T a, T b) {
    a + b; a += b;
    a - b; a -= b;
    a * b; a *= b;
    a / b; a /= b;
};

template<bool enable>
struct invokeIf {};

template<>
struct invokeIf<true>
{
    template<Invokable FunctionType>
    static inline void invoke(FunctionType function)
    { function(); }
};

template<>
struct invokeIf<false>
{
    template<Invokable FunctionType>
    static inline void invoke(FunctionType function) {}
};