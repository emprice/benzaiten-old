#ifndef _BZEXPRESSION_HH_
#define _BZEXPRESSION_HH_

#include <vector>
#include <string>
#include <unordered_map>

namespace benzaiten
{
    template <typename E>
    struct FunctionExpression
    {
    };

    struct SubstituteEntry
    {
        SubstituteEntry(const std::string &name, double value,
            const std::unordered_map<std::string, size_t> &d) :
                name(name), value(value), d(d)
        {
            // no-op
        }

        ~SubstituteEntry() { }

        std::string name;
        double value;
        std::unordered_map<std::string, size_t> d;
    };
}

#endif      // _BZEXPRESSION_HH_

// vim: set ft=cpp.doxygen:
