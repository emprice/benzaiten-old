#ifndef _BZVARIABLE_HH_
#define _BZVARIABLE_HH_

#include "bzexpression.hh"

#include <string>
#include <iostream>

namespace benzaiten
{
    enum VariableType
    {
        Spatial,
        Temporal
    };

    struct Variable : public FunctionExpression<Variable>
    {
        template <size_t Order>
        using deriv_type = Variable;

        Variable(const std::string &name, VariableType type) :
            name(name), type(type) { }

        Variable(const Variable &other) : name(other.name),
            type(other.type), _isConcrete(other._isConcrete),
                _value(other._value) { }

        std::string getName() const { return name; }

        VariableType getType() const { return type; }

        template <size_t Order = 1>
        Variable& derivativeInPlace(const Variable &other)
        {
            if (name == other.name)
            {
                _isConcrete = true;
                _value = 1;
            }
            else
            {
                _isConcrete = true;
                _value = 0;
            }

            return *this;
        }

        template <size_t Order = 1>
        Variable derivative(const Variable &other) const
        {
            return Variable(*this).derivativeInPlace<Order>(other);
        }

        bool isConcrete() const { return _isConcrete; }

        double getValue() const { return _value; }

        friend std::ostream& operator<<(std::ostream &os, const Variable &vbl)
        {
            if (vbl._isConcrete) os << vbl._value;
            else os << vbl.name;

            return os;
        }

        private:
            std::string name;
            VariableType type;

            bool _isConcrete = false;
            double _value;
    };
}

#endif      // _BZVARIABLE_HH_

// vim: set ft=cpp.doxygen:
