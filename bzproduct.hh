#ifndef _BZPRODUCT_HH_
#define _BZPRODUCT_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzsum.hh"

namespace benzaiten
{
    template <typename E1, typename E2>
    struct FunctionProduct;

    template <typename F1, typename F2, size_t Order>
    struct ProductDerivativeType
    {
        using type =
            FunctionSum<typename ProductDerivativeType<typename F1::template deriv_type<1>, F2, Order - 1>::type,
                typename ProductDerivativeType<F1, typename F2::template deriv_type<1>, Order - 1>::type>;
    };

    template <typename F1, typename F2>
    struct ProductDerivativeType<F1, F2, 0>
    {
        using type = FunctionProduct<F1, F2>;
    };

    template <typename E1, typename E2>
    struct FunctionProduct : public FunctionExpression<FunctionProduct<E1, E2>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<E1, E2, Order>::type;

            FunctionProduct(const E1 &_fn1, const E2 &_fn2) : fn1(_fn1), fn2(_fn2) { }

            template <size_t Order = 1>
            typename ProductDerivativeType<E1, E2, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (fn1.template derivative<1>(var) * fn2).template derivative<Order-1>(var) +
                            (fn1 * fn2.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionProduct<E1, E2>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn1.substituteInPlace(subs);
                fn2.substituteInPlace(subs);

                if (fn1.isConcrete() && fn2.isConcrete())
                {
                    _isConcrete = true;
                    _value = fn1.getValue() * fn2.getValue();
                }

                return *this;
            }

            FunctionProduct<E1, E2> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionProduct<E1, E2>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionProduct<E1, E2> &prod)
            {
                if (prod._isConcrete) os << prod._value;
                else os << "(" << prod.fn1 << " * " << prod.fn2 << ")";

                return os;
            }

        private:
            E1 fn1;
            E2 fn2;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E1, typename E2>
    FunctionProduct<E1, E2> operator*(FunctionExpression<E1> const& fn1,
        FunctionExpression<E2> const& fn2)
    {
        return FunctionProduct<E1, E2>(static_cast<E1 const&>(fn1),
            static_cast<E2 const&>(fn2));
    }

    template <typename E>
    FunctionProduct<E, Constant> operator*(FunctionExpression<E> const& fn, double other)
    {
        return FunctionProduct<E, Constant>(static_cast<E const&>(fn), Constant(other));
    }

    template <typename E>
    FunctionProduct<Constant, E> operator*(double other, FunctionExpression<E> const& fn)
    {
        return FunctionProduct<Constant, E>(Constant(other), static_cast<E const&>(fn));
    }
}

#endif      // _BZPRODUCT_HH_

// vim: set ft=cpp.doxygen:
