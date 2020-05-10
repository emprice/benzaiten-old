#ifndef _BZQUOTIENT_HH_
#define _BZQUOTIENT_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzsum.hh"

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionQuotient;

    template <typename F1, typename F2, size_t Order>
    struct QuotientDerivativeType
    {
        using type = FunctionDifference<FunctionQuotient<typename F1::template deriv_type<1>, F2>,
            FunctionQuotient<FunctionProduct<F1, typename F2::template deriv_type<1>>, FunctionProduct<F2, F2>>>;
    };

    template <typename F1, typename F2>
    struct QuotientDerivativeType<F1, F2, 0>
    {
        using type = FunctionQuotient<F1, F2>;
    };

    template <typename E1, typename E2>
    struct FunctionQuotient : public FunctionExpression<FunctionQuotient<E1, E2>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename QuotientDerivativeType<E1, E2, Order>::type;

            FunctionQuotient(const E1 &fn1, const E2 &fn2) : fn1(fn1), fn2(fn2) { }

            template <size_t Order = 1>
            typename QuotientDerivativeType<E1, E2, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return ((fn1.template derivative<1>(var) / fn2) -
                    (fn1 * fn2.template derivative<1>(var) / (fn2 * fn2))).template derivative<Order-1>(var);
            }

            FunctionQuotient<E1, E2>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);
                fn2.substituteInPlace(subs);

                if (fn1.isConcrete() && fn2.isConcrete())
                {
                    _isConcrete = true;
                    _value = fn1.getValue() / fn2.getValue();
                }

                return *this;
            }

            FunctionQuotient<E1, E2> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionQuotient<E1, E2>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionQuotient<E1, E2> &quo)
            {
                if (quo._isConcrete) os << quo._value;
                else os << "(" << quo.fn1 << " / " << quo.fn2 << ")";

                return os;
            }

        private:
            E1 fn1;
            E2 fn2;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E1, typename E2>
    FunctionQuotient<E1, E2> operator/(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionQuotient<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }
}

#endif      // _BZQUOTIENT_HH_

// vim: set ft=cpp.doxygen:
