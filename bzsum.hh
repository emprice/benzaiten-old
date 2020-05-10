#ifndef _BZSUM_HH_
#define _BZSUM_HH_

#include "bzfunction.hh"

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionSum;

    template <typename F1, typename F2, size_t Order>
    struct SumDerivativeType
    {
        using type = FunctionSum<typename F1::template deriv_type<Order>,
            typename F2::template deriv_type<Order>>;
    };

    template <typename F1, typename F2>
    struct SumDerivativeType<F1, F2, 0>
    {
        using type = FunctionSum<F1, F2>;
    };

    template <typename E1, typename E2>
    struct FunctionSum : public FunctionExpression<FunctionSum<E1, E2>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = FunctionSum<typename E1::template deriv_type<Order>,
                typename E2::template deriv_type<Order>>;

            FunctionSum(const E1 &fn1, const E2 &fn2) : fn1(fn1), fn2(fn2) { }

            template <size_t Order = 1>
            typename SumDerivativeType<E1, E2, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return fn1.template derivative<Order>(var) + fn2.template derivative<Order>(var);
            }

            FunctionSum<E1, E2>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);
                fn2.substituteInPlace(subs);

                if (fn1.isConcrete() && fn2.isConcrete())
                {
                    _isConcrete = true;
                    _value = fn1.getValue() + fn2.getValue();
                }

                return *this;
            }

            FunctionSum<E1, E2> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionSum<E1, E2>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionSum<E1, E2> &sum)
            {
                if (sum._isConcrete) os << sum._value;
                else os << "(" << sum.fn1 << " + " << sum.fn2 << ")";

                return os;
            }

        private:
            E1 fn1;
            E2 fn2;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E1, typename E2>
    FunctionSum<E1, E2> operator+(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionSum<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E>
    FunctionSum<E, Constant> operator+(FunctionExpression<E> const& fn, double other)
    {
        return FunctionSum<E, Constant>(static_cast<E const&>(fn), Constant(other));
    }

    template <typename E>
    FunctionSum<Constant, E> operator+(double other, FunctionExpression<E> const& fn)
    {
        return FunctionSum<Constant, E>(Constant(other), static_cast<E const&>(fn));
    }
}

#endif      // _BZSUM_HH_

// vim: set ft=cpp.doxygen:
