#ifndef _BZPOWER_HH_
#define _BZPOWER_HH_

#include <cmath>

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
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionPower<E1, FunctionDifference<E2, Constant>>,
                FunctionSum<FunctionProduct<E2, typename E1::template deriv_type<1>>,
                    FunctionProduct<FunctionProduct<E1, typename E2::template deriv_type<1>>,
                FunctionLog<E1>>>, Order - 1>::type;

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

            FunctionPower<E1, E2>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);
                fn2.substituteInPlace(subs);

                if (fn1.isConcrete() && fn2.isConcrete())
                {
                    _isConcrete = true;
                    _value = pow(fn1.getValue(), fn2.getValue());
                }

                return *this;
            }

            FunctionPower<E1, E2> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionPower<E1, E2>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionPower<E1, E2> &pwr)
            {
                if (pwr._isConcrete) os << pwr._value;
                else os << "(" << pwr.fn1 << " ^ " << pwr.fn2 << ")";

                return os;
            }

        private:
            E1 fn1;
            E2 fn2;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E1, typename E2>
    FunctionPower<E1, E2> operator^(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionPower<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E>
    FunctionPower<E, Constant> operator^(FunctionExpression<E> const& fn, double pwr)
    {
        return FunctionPower<E, Constant>(static_cast<E const&>(fn), Constant(pwr));
    }

    template <typename E1, typename E2>
    FunctionPower<E1, E2> pow(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionPower<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E>
    FunctionPower<E, Constant> pow(FunctionExpression<E> const& fn, double pwr)
    {
        return FunctionPower<E, Constant>(static_cast<E const&>(fn), Constant(pwr));
    }

    template <typename E>
    FunctionPower<E, Constant> sqrt(FunctionExpression<E> const& fn)
    {
        return FunctionPower<E, Constant>(static_cast<E const&>(fn), Constant(0.5));
    }
}

#endif      // _BZPOWER_HH_

// vim: set ft=cpp.doxygen:
