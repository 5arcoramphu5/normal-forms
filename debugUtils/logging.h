#ifndef _LOGGING_
#define _LOGGING_

#include <string>

enum VerbosityLevel{None, Minimal, Diagnostic};

template<typename T>
concept LoggerType = requires (std::string message) {
    { T::template print<VerbosityLevel::None>(message) };
};

template <typename T>
concept Streamable = requires(std::ostream &os, T value) {
    { os << value } -> std::convertible_to<std::ostream &>;
};

template<VerbosityLevel Verbosity>
class Logger
{
    public:
    template<VerbosityLevel MessageVerbosity, Streamable MessageType>
    static void print(MessageType message)
    {
        static_assert(MessageVerbosity != VerbosityLevel::None);

        if constexpr (Verbosity >= MessageVerbosity)
            std::cout << message << std::endl;
    }
};

#endif