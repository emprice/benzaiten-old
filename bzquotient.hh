#ifndef _BZQUOTIENT_HH_
#define _BZQUOTIENT_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzneg.hh"
#include "bzproduct.hh"

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionQuotient;

    template <typename E1>
    struct FunctionQuotientSimple1;

    template <typename E2>
    struct FunctionQuotientSimple2;

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

    template <typename F1, size_t Order>
    struct Simple1QuotientDerivativeType
    {
        using type = FunctionQuotientSimple1<typename F1::template deriv_type<Order>>;
    };

    template <typename F1>
    struct Simple1QuotientDerivativeType<F1, 0>
    {
        using type = FunctionQuotientSimple1<F1>;
    };

    template <typename F2, size_t Order>
    struct Simple2QuotientDerivativeType
    {
        using type = FunctionProductSimple<typename QuotientDerivativeType<FunctionNegate<typename F2::template deriv_type<1>>, FunctionProduct<F2, F2>, Order - 1>::type>;
    };

    template <typename F2>
    struct Simple2QuotientDerivativeType<F2, 0>
    {
        using type = FunctionQuotientSimple2<F2>;
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

    template <typename E1>
    struct FunctionQuotientSimple1 : public FunctionExpression<FunctionQuotientSimple1<E1>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename Simple1QuotientDerivativeType<E1, Order>::type;

            FunctionQuotientSimple1(const E1 &fn1, const Constant &cnst) : fn1(fn1), cnst(cnst) { }

            template <size_t Order = 1>
            typename Simple1QuotientDerivativeType<E1, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return fn1.template derivative<Order>(var) / cnst.getValue();
            }

            FunctionQuotientSimple1<E1>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);

                if (fn1.isConcrete())
                {
                    _isConcrete = true;
                    _value = fn1.getValue() / cnst.getValue();
                }

                return *this;
            }

            FunctionQuotientSimple1<E1> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionQuotientSimple1<E1>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionQuotientSimple1<E1> &quo)
            {
                if (quo._isConcrete) os << quo._value;
                else os << "(" << quo.fn1 << " / " << quo.cnst << ")";

                return os;
            }

        private:
            E1 fn1;
            Constant cnst;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionQuotientSimple1<E> operator/(FunctionExpression<E> const& fn, double other)
    {
        return FunctionQuotientSimple1<E>(static_cast<E const&>(fn), Constant(other));
    }

    template <typename E2>
    struct FunctionQuotientSimple2 : public FunctionExpression<FunctionQuotientSimple2<E2>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename Simple2QuotientDerivativeType<E2, Order>::type;

            FunctionQuotientSimple2(const Constant &cnst, const E2 &fn2) : cnst(cnst), fn2(fn2) { }

            template <size_t Order = 1>
            typename Simple2QuotientDerivativeType<E2, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return  cnst.getValue() * (-fn2.template derivative<1>(var) /
                    (fn2 * fn2)).template derivative<Order - 1>(var);
            }

            FunctionQuotientSimple2<E2>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                if (cnst.getValue() == 0)
                {
                    _isConcrete = true;
                    _value = 0;
                }
                else
                {
                    fn2.substituteInPlace(subs);

                    if (fn2.isConcrete())
                    {
                        _isConcrete = true;
                        _value = cnst.getValue() / fn2.getValue();
                    }
                }

                return *this;
            }

            FunctionQuotientSimple2<E2> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionQuotientSimple2<E2>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionQuotientSimple2<E2> &quo)
            {
                if (quo._isConcrete) os << quo._value;
                else os << "(" << quo.cnst << " / " << quo.fn2 << ")";

                return os;
            }

        private:
            E2 fn2;
            Constant cnst;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionQuotientSimple2<E> operator/(double other, FunctionExpression<E> const& fn)
    {
        return FunctionQuotientSimple2<E>(Constant(other), static_cast<E const&>(fn));
    }
}

#endif      // _BZQUOTIENT_HH_

// vim: set ft=cpp.doxygen:
