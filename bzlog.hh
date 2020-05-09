#ifndef _BZLOG_HH_
#define _BZLOG_HH_

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
        using type = QuotientDerivativeType<F, typename F::template deriv_type<1>, Order - 1>;
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
            FunctionLog(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename LogDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (fn.template derivative<1>(var) / fn).template derivative<Order-1>(var);
            }

            friend std::ostream& operator<<(std::ostream &os, const FunctionLog<E> &expr)
            {
                os << "log(" << expr.fn << ")";
                return os;
            }

        private:
            const E fn;
    };

    template <typename E>
    FunctionLog<E> log(FunctionExpression<E> const& fn)
    {
        return FunctionLog<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZLOG_HH_

// vim: set ft=cpp.doxygen: