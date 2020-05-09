#ifndef _BZPOWER_HH_
#define _BZPOWER_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzproduct.hh"
#include "bzlog.hh"

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionPower;

    template <typename F1, typename F2, size_t Order>
    struct PowerDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionPower<F1, FunctionDifference<F2, Constant>>,
            FunctionSum<FunctionProduct<F2, typename F1::template deriv_type<1>>,
            FunctionProduct<FunctionProduct<F1, typename F2::template deriv_type<1>>,
                FunctionLog<F1>>>, Order - 1>::type;
    };

    template <typename F1, typename F2>
    struct PowerDerivativeType<F1, F2, 0>
    {
        using type = FunctionPower<F1, F2>;
    };

    template <typename E1, typename E2>
    struct FunctionPower : public FunctionExpression<FunctionPower<E1, E2>>
    {
        public:
            FunctionPower(const E1 &fn1, const E2 &fn2) : fn1(fn1), fn2(fn2) { }

            template <size_t Order = 1>
            typename PowerDerivativeType<E1, E2, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else
                {
                    Constant one(1);
                    return (pow(fn1, fn2 - one) * (fn2 * fn1.template derivative<1>(var) +
                        fn1 * fn2.template derivative<1>(var) * log(fn1))).template derivative<Order-1>(var);
                }
            }

            friend std::ostream& operator<<(std::ostream &os, const FunctionPower<E1, E2> &pwr)
            {
                os << "(" << pwr.fn1 << " ^ " << pwr.fn2 << ")";
                return os;
            }

        private:
            const E1 fn1;
            const E2 fn2;
    };

    template <typename E1, typename E2>
    FunctionPower<E1, E2> operator^(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionPower<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E1, typename E2>
    FunctionPower<E1, E2> pow(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionPower<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }
}

#endif      // _BZPOWER_HH_

// vim: set ft=cpp.doxygen:
