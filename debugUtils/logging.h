#ifndef _LOGGING_
#define _LOGGING_

#include <string>

enum VerbosityLevel{None, Minimal, Diagnostic};

template<VerbosityLevel Verbosity>
class Logger
{
    public:
    template<VerbosityLevel MessageVerbosity>
    static void print(std::string message)
    {
        static_assert(MessageVerbosity != VerbosityLevel::None);

        if constexpr (Verbosity >= MessageVerbosity)
            std::cout << message << std::endl;
    }
};

template<typename T>
concept LoggerType = requires (std::string message) {
    { T::template print<VerbosityLevel::None>(message) };
};

#endif