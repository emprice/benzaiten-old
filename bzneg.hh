#ifndef _BZNEG_HH_
#define _BZNEG_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"

namespace benzaiten
{
    template <typename E>
    struct FunctionNegate;

    template <typename E>
    struct FunctionNegate : public FunctionExpression<FunctionNegate<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = FunctionNegate<typename E::template deriv_type<Order>>;

            FunctionNegate(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            FunctionNegate<typename E::template deriv_type<Order>>
                derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return -(fn.template derivative<Order>(var));
            }

            FunctionNegate<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = -(fn.getValue());
                }

                return *this;
            }

            FunctionNegate<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionNegate<E>(*this).substitute(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionNegate<E> &neg)
            {
                if (neg._isConcrete) os << neg._value;
                else os << "(-" << neg.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionNegate<E> operator-(FunctionExpression<E> const& fn)
    {
        return FunctionNegate<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZNEG_HH_

// vim: set ft=cpp.doxygen:
