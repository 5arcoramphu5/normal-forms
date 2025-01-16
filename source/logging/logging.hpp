#pragma once

#include <string>
#include "capd/capdlib.h"
#include "../templateUtils.hpp"
#include "polynomialPrintingPolicies.hpp"

enum VerbosityLevel{None, Minimal, Diagnostic, Debug, Error};

template<typename T>
concept LoggerType = requires (std::string message) {
    { T::template log<Minimal>(message) };
};

template<VerbosityLevel Verbosity, PolynomialPrintingPolicy PolynomialPrinting = SymbolicPolynomialPrinting>
class Logger
{
    public:

    template<VerbosityLevel MessageVerbosity>
    static void log() 
    { 
        invokeIf<Verbosity >= MessageVerbosity>::template invoke([]()
        { std::cout << std::endl; });
    }

    template<VerbosityLevel MessageVerbosity, Streamable MessageType, Streamable... MessageTypes>
    static void log(const MessageType &message, MessageTypes... messages)
    {
        static_assert(MessageVerbosity != None, "cannot log message with verbosity None");

        invokeIf<Verbosity >= MessageVerbosity>::template invoke([message, messages...]()
        {
            std::cout << message << " ";
            log<MessageVerbosity>(messages...);
        });
    }

    template<VerbosityLevel MessageVerbosity, ArithmeticType Coeff, Streamable... MessageTypes>
    static void log(const Polynomial<Coeff> &polynomial, MessageTypes... messages)
    {
        static_assert(MessageVerbosity != None, "cannot log message with verbosity None");

        invokeIf<Verbosity >= MessageVerbosity>::template invoke([polynomial, messages...]()
        {
            std::cout << PolynomialPrinting::polyToString(polynomial) << " ";
            log<MessageVerbosity>(messages...);
        });
    }

    // enable function call based on verbosity
    template<VerbosityLevel EnablingVerbosity, Invokable FunctionType>
    static inline void enableIf(const FunctionType &function)
    { 
        invokeIf<Verbosity >= EnablingVerbosity>::template invoke(function);
    }

};