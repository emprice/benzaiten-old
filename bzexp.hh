#ifndef _BZEXP_HH_
#define _BZEXP_HH_

#include <cmath>

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzproduct.hh"

namespace benzaiten
{
    template <typename E>
    struct FunctionExp;

    template <typename F, size_t Order>
    struct ExpDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionExp<F>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct ExpDerivativeType<F, 0>
    {
        using type = FunctionExp<F>;
    };

    template <typename E>
    struct FunctionExp : public FunctionExpression<FunctionExp<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionExp<E>,
                typename E::template deriv_type<1>, Order - 1>::type;

            FunctionExp(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename ExpDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (exp(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionExp<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = exp(fn.getValue());
                }

                return *this;
            }

            FunctionExp<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionExp<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionExp<E> &expr)
            {
                os << "exp(" << expr.fn << ")";
                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionExp<E> exp(FunctionExpression<E> const& fn)
    {
        return FunctionExp<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZEXP_HH_

// vim: set ft=cpp.doxygen:
