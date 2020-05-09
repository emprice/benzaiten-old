#ifndef _BZVARIABLE_HH_
#define _BZVARIABLE_HH_

#include <string>

namespace benzaiten
{
    enum VariableType
    {
        Spatial,
        Temporal
    };

    struct Variable
    {
        Variable(const std::string &name, VariableType type) :
            name(name), type(type) { }

        std::string getName() const { return name; }

        VariableType getType() const { return type; }

        private:
            std::string name;
            VariableType type;
    };
}

#endif      // _BZVARIABLE_HH_

// vim: set ft=cpp.doxygen:
