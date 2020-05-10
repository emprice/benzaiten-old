#ifndef _BZLOG_HH_
#define _BZLOG_HH_

#include <cmath>

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzquotient.hh"

namespace benzaiten
{
    template <typename E>
    struct FunctionLog;

    template <typename F, size_t Order>
    struct LogDerivativeType
    {
        using type = typename QuotientDerivativeType<F, typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct LogDerivativeType<F, 0>
    {
        using type = FunctionLog<F>;
    };

    template <typename E>
    struct FunctionLog : public FunctionExpression<FunctionLog<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename QuotientDerivativeType<E, typename E::template deriv_type<1>, Order - 1>::type;

            FunctionLog(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename LogDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (fn.template derivative<1>(var) / fn).template derivative<Order-1>(var);
            }

            FunctionLog<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = log(fn.getValue());
                }

                return *this;
            }

            FunctionLog<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionLog<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionLog<E> &expr)
            {
                os << "log(" << expr.fn << ")";
                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionLog<E> log(FunctionExpression<E> const& fn)
    {
        return FunctionLog<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZLOG_HH_

// vim: set ft=cpp.doxygen:
