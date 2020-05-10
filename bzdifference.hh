#ifndef _BZDIFFERENCE_HH_
#define _BZDIFFERENCE_HH_

#include "bzfunction.hh"

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionDifference;

    template <typename F1, typename F2, size_t Order>
    struct DifferenceDerivativeType
    {
        using type = FunctionDifference<typename F1::template deriv_type<Order>,
            typename F2::template deriv_type<Order>>;
    };

    template <typename F1, typename F2>
    struct DifferenceDerivativeType<F1, F2, 0>
    {
        using type = FunctionDifference<F1, F2>;
    };

    template <typename E1, typename E2>
    struct FunctionDifference : public FunctionExpression<FunctionDifference<E1, E2>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = FunctionDifference<typename E1::template deriv_type<Order>,
                  typename E2::template deriv_type<Order>>;

            FunctionDifference(const E1 &fn1, const E2 &fn2) : fn1(fn1), fn2(fn2) { }

            template <size_t Order = 1>
            typename DifferenceDerivativeType<E1, E2, Order>::type derivative(const Variable &var) const
            {
                return fn1.template derivative<Order>(var) - fn2.template derivative<Order>(var);
            }

            FunctionDifference<E1, E2>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);
                fn2.substituteInPlace(subs);

                if (fn1.isConcrete() && fn2.isConcrete())
                {
                    _isConcrete = true;
                    _value = fn1.getValue() - fn2.getValue();
                }

                return *this;
            }

            FunctionDifference<E1, E2> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionDifference<E1, E2>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionDifference<E1, E2> &diff)
            {
                if (diff._isConcrete) os << diff._value;
                else os << "(" << diff.fn1 << " - " << diff.fn2 << ")";

                return os;
            }

        private:
            E1 fn1;
            E2 fn2;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E1, typename E2>
    FunctionDifference<E1, E2> operator-(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionDifference<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E>
    FunctionDifference<E, Constant> operator-(FunctionExpression<E> const& fn, double other)
    {
        return FunctionDifference<E, Constant>(static_cast<E const&>(fn), Constant(other));
    }

    template <typename E>
    FunctionDifference<Constant, E> operator-(double other, FunctionExpression<E> const& fn)
    {
        return FunctionDifference<Constant, E>(Constant(other), static_cast<E const&>(fn));
    }
}

#endif      // _BZDIFFERENCE_HH_

// vim: set ft=cpp.doxygen:
