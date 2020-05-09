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

            friend std::ostream& operator<<(std::ostream &os, const FunctionSum<E1, E2> &sum)
            {
                os << "(" << sum.fn1 << " + " << sum.fn2 << ")";
                return os;
            }

        private:
            const E1 fn1;
            const E2 fn2;
    };

    template <typename E1, typename E2>
    FunctionSum<E1, E2> operator+(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionSum<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }
}

#endif      // _BZSUM_HH_

// vim: set ft=cpp.doxygen: