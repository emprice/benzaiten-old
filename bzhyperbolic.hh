#ifndef _BZHYPERBOLIC_HH_
#define _BZHYPERBOLIC_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"

#include <cmath>

double coth(double x)
{
    return cosh(x) / sinh(x);
}

double sech(double x)
{
    return 1. / cosh(x);
}

double csch(double x)
{
    return 1. / sinh(x);
}

namespace benzaiten
{
    template <typename E>
    struct FunctionSinh;

    template <typename E>
    struct FunctionCosh;

    template <typename E>
    struct FunctionTanh;

    template <typename E>
    struct FunctionCoth;

    template <typename E>
    struct FunctionSech;

    template <typename E>
    struct FunctionCsch;

    template <typename F, size_t Order>
    struct SinhDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionCosh<F>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct SinhDerivativeType<F, 0>
    {
        using type = FunctionSinh<F>;
    };

    template <typename F, size_t Order>
    struct CoshDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionSinh<F>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct CoshDerivativeType<F, 0>
    {
        using type = FunctionCosh<F>;
    };

    template <typename F, size_t Order>
    struct TanhDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionSech<F>, FunctionSech<F>>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct TanhDerivativeType<F, 0>
    {
        using type = FunctionTanh<F>;
    };

    template <typename F, size_t Order>
    struct CothDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCsch<F>>, FunctionCsch<F>>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct CothDerivativeType<F, 0>
    {
        using type = FunctionCoth<F>;
    };

    template <typename F, size_t Order>
    struct SechDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionTanh<F>>, FunctionSech<F>>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct SechDerivativeType<F, 0>
    {
        using type = FunctionSech<F>;
    };

    template <typename F, size_t Order>
    struct CschDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCoth<F>>, FunctionCsch<F>>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct CschDerivativeType<F, 0>
    {
        using type = FunctionCsch<F>;
    };

    template <typename E>
    struct FunctionSinh : public FunctionExpression<FunctionSinh<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionCosh<E>,
                typename E::template deriv_type<1>, Order - 1>::type;

            FunctionSinh(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename SinhDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (cosh(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionSinh<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = sinh(fn.getValue());
                }

                return *this;
            }

            FunctionSinh<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionSinh<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionSinh<E> &hyper)
            {
                if (hyper._isConcrete) os << hyper._value;
                else os << "sinh(" << hyper.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionSinh<E> sinh(FunctionExpression<E> const& fn)
    {
        return FunctionSinh<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionCosh : public FunctionExpression<FunctionCosh<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionSinh<E>,
                typename E::template deriv_type<1>, Order - 1>::type;

            FunctionCosh(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename CoshDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (sinh(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionCosh<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = sinh(fn.getValue());
                }

                return *this;
            }

            FunctionCosh<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionCosh<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionCosh<E> &hyper)
            {
                if (hyper._isConcrete) os << hyper._value;
                else os << "cosh(" << hyper.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionCosh<E> cosh(FunctionExpression<E> const& fn)
    {
        return FunctionCosh<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionTanh : public FunctionExpression<FunctionTanh<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionProduct<FunctionSech<E>, FunctionSech<E>>,
                typename E::template deriv_type<1>, Order - 1>::type;

            FunctionTanh(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename TanhDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (sech(fn) * sech(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionTanh<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = tanh(fn.getValue());
                }

                return *this;
            }

            FunctionTanh<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionTanh<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionTanh<E> &hyper)
            {
                if (hyper._isConcrete) os << hyper._value;
                else os << "tanh(" << hyper.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionTanh<E> tanh(FunctionExpression<E> const& fn)
    {
        return FunctionTanh<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionCoth : public FunctionExpression<FunctionCoth<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type =
                typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCsch<E>>, FunctionCsch<E>>,
                    typename E::template deriv_type<1>, Order - 1>::type;

            FunctionCoth(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename CothDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (-csch(fn) * csch(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionCoth<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = coth(fn.getValue());
                }

                return *this;
            }

            FunctionCoth<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionCoth<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionCoth<E> &hyper)
            {
                if (hyper._isConcrete) os << hyper._value;
                else os << "coth(" << hyper.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionCoth<E> coth(FunctionExpression<E> const& fn)
    {
        return FunctionCoth<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionSech : public FunctionExpression<FunctionSech<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type =
                typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionTanh<E>>, FunctionSech<E>>,
                    typename E::template deriv_type<1>, Order - 1>::type;

            FunctionSech(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename SechDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (-tanh(fn) * sech(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionSech<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = sech(fn.getValue());
                }

                return *this;
            }

            FunctionSech<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionSech<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionSech<E> &hyper)
            {
                if (hyper._isConcrete) os << hyper._value;
                else os << "sech(" << hyper.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionSech<E> sech(FunctionExpression<E> const& fn)
    {
        return FunctionSech<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionCsch : public FunctionExpression<FunctionCsch<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type =
                typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCoth<E>>, FunctionCsch<E>>,
                    typename E::template deriv_type<1>, Order - 1>::type;

            FunctionCsch(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename CschDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (-coth(fn) * csch(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionCsch<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = csch(fn.getValue());
                }

                return *this;
            }

            FunctionCsch<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionCsch<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionCsch<E> &hyper)
            {
                if (hyper._isConcrete) os << hyper._value;
                else os << "csch(" << hyper.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionCsch<E> csch(FunctionExpression<E> const& fn)
    {
        return FunctionCsch<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZHYPERBOLIC_HH_

// vim: set ft=cpp.doxygen:
