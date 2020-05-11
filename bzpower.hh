#ifndef _BZPOWER_HH_
#define _BZPOWER_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzproduct.hh"
#include "bzlog.hh"

#include <cmath>

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionPower;

    template <typename E1>
    struct FunctionPowerSimple;

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

    template <typename F, size_t Order>
    struct SimplePowerDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<Constant,
            FunctionPowerSimple<F>>, typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct SimplePowerDerivativeType<F, 0>
    {
        using type = FunctionPowerSimple<F>;
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

    template <typename E1, typename E2>
    FunctionPower<E1, E2> pow(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionPower<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E1>
    struct FunctionPowerSimple : public FunctionExpression<FunctionPowerSimple<E1>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionProduct<Constant, FunctionPowerSimple<E1>>,
                typename E1::template deriv_type<1>, Order - 1>::type;

            FunctionPowerSimple(const E1 &fn1, const Constant &cnst) : fn1(fn1), cnst(cnst) { }

            template <size_t Order = 1>
            typename SimplePowerDerivativeType<E1, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else
                {
                    return (cnst * pow(fn1, cnst.getValue() - 1) *
                        fn1.template derivative<1>(var)).template derivative<Order-1>(var);
                }
            }

            FunctionPowerSimple<E1>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);

                if (fn1.isConcrete())
                {
                    _isConcrete = true;
                    _value = std::pow(fn1.getValue(), cnst.getValue());
                }

                return *this;
            }

            FunctionPowerSimple<E1> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionPowerSimple<E1>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionPowerSimple<E1> &pwr)
            {
                if (pwr._isConcrete) os << pwr._value;
                else os << "(" << pwr.fn1 << " ^ " << pwr.cnst.getValue() << ")";

                return os;
            }

        private:
            E1 fn1;
            Constant cnst;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionPowerSimple<E> operator^(FunctionExpression<E> const& fn, double pwr)
    {
        return FunctionPowerSimple<E>(static_cast<E const&>(fn), Constant(pwr));
    }

    template <typename E>
    FunctionPowerSimple<E> pow(FunctionExpression<E> const& fn, double pwr)
    {
        return FunctionPowerSimple<E>(static_cast<E const&>(fn), Constant(pwr));
    }

    template <typename E>
    FunctionPowerSimple<E> sqrt(FunctionExpression<E> const& fn)
    {
        return FunctionPowerSimple<E>(static_cast<E const&>(fn), Constant(0.5));
    }
}

#endif      // _BZPOWER_HH_

// vim: set ft=cpp.doxygen:
