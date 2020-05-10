#ifndef _BZTRIG_HH_
#define _BZTRIG_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzneg.hh"

#include <cmath>

double sec(double x)
{
    return 1. / cos(x);
}

double csc(double x)
{
    return 1. / sin(x);
}

double cot(double x)
{
    return 1. / tan(x);
}

namespace benzaiten
{
    template <typename E>
    struct FunctionSine;

    template <typename E>
    struct FunctionCosine;

    template <typename E>
    struct FunctionTangent;

    template <typename E>
    struct FunctionCotangent;

    template <typename E>
    struct FunctionSecant;

    template <typename E>
    struct FunctionCosecant;

    template <typename F, size_t Order>
    struct SineDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionCosine<F>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct SineDerivativeType<F, 0>
    {
        using type = FunctionSine<F>;
    };

    template <typename F, size_t Order>
    struct CosineDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionNegate<FunctionSine<F>>,
            typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct CosineDerivativeType<F, 0>
    {
        using type = FunctionCosine<F>;
    };

    template <typename F, size_t Order>
    struct TangentDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionSecant<F>,
            FunctionSecant<F>>, typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct TangentDerivativeType<F, 0>
    {
        using type = FunctionTangent<F>;
    };

    template <typename F, size_t Order>
    struct CotangentDerivativeType
    {
        using type =
            typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCosecant<F>>,
                FunctionCosecant<F>>, typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct CotangentDerivativeType<F, 0>
    {
        using type = FunctionCotangent<F>;
    };

    template <typename F, size_t Order>
    struct SecantDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionSecant<F>,
            FunctionTangent<F>>, typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct SecantDerivativeType<F, 0>
    {
        using type = FunctionSecant<F>;
    };

    template <typename F, size_t Order>
    struct CosecantDerivativeType
    {
        using type = typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCosecant<F>>,
            FunctionCotangent<F>>, typename F::template deriv_type<1>, Order - 1>::type;
    };

    template <typename F>
    struct CosecantDerivativeType<F, 0>
    {
        using type = FunctionCosecant<F>;
    };

    template <typename E>
    struct FunctionSine : public FunctionExpression<FunctionSine<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionCosine<E>,
                typename E::template deriv_type<1>, Order - 1>::type;

            FunctionSine(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename SineDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (cos(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionSine<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = sin(fn.getValue());
                }

                return *this;
            }

            FunctionSine<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionSine<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionSine<E> &sine)
            {
                if (sine._isConcrete) os << sine._value;
                else os << "sin(" << sine.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionSine<E> sin(FunctionExpression<E> const& fn)
    {
        return FunctionSine<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionCosine : public FunctionExpression<FunctionCosine<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionNegate<FunctionSine<E>>,
                typename E::template deriv_type<1>, Order - 1>::type;

            FunctionCosine(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename CosineDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (-sin(fn) * fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionCosine<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = cos(fn.getValue());
                }

                return *this;
            }

            FunctionCosine<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionCosine<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionCosine<E> &cosine)
            {
                if (cosine._isConcrete) os << cosine._value;
                else os << "cos(" << cosine.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionCosine<E> cos(FunctionExpression<E> const& fn)
    {
        return FunctionCosine<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionTangent : public FunctionExpression<FunctionTangent<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionProduct<FunctionSecant<E>,
                FunctionSecant<E>>, typename E::template deriv_type<1>, Order - 1>::type;

            FunctionTangent(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename TangentDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (sec(fn) * sec(fn) *
                    fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionTangent<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = tan(fn.getValue());
                }

                return *this;
            }

            FunctionTangent<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionTangent<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionTangent<E> &tngt)
            {
                if (tngt._isConcrete) os << tngt._value;
                else os << "tan(" << tngt.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionTangent<E> tan(FunctionExpression<E> const& fn)
    {
        return FunctionTangent<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionCotangent : public FunctionExpression<FunctionCotangent<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCosecant<E>>,
                FunctionCosecant<E>>, typename E::template deriv_type<1>, Order - 1>::type;

            FunctionCotangent(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename CotangentDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (-csc(fn) * csc(fn) *
                    fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionCotangent<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = cot(fn.getValue());
                }

                return *this;
            }

            FunctionCotangent<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionCotangent<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionCotangent<E> &cotan)
            {
                if (cotan._isConcrete) os << cotan._value;
                else os << "cot(" << cotan.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionCotangent<E> cot(FunctionExpression<E> const& fn)
    {
        return FunctionCotangent<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionSecant : public FunctionExpression<FunctionSecant<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionProduct<FunctionSecant<E>,
                FunctionTangent<E>>, typename E::template deriv_type<1>, Order - 1>::type;

            FunctionSecant(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename SecantDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (sec(fn) * tan(fn) *
                    fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionSecant<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = sec(fn.getValue());
                }

                return *this;
            }

            FunctionSecant<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionSecant<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionSecant<E> &secant)
            {
                if (secant._isConcrete) os << secant._value;
                else os << "sec(" << secant.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionSecant<E> sec(FunctionExpression<E> const& fn)
    {
        return FunctionSecant<E>(static_cast<E const&>(fn));
    }

    template <typename E>
    struct FunctionCosecant : public FunctionExpression<FunctionCosecant<E>>
    {
        public:
            template <size_t Order = 1>
            using deriv_type = typename ProductDerivativeType<FunctionProduct<FunctionNegate<FunctionCosecant<E>>,
                FunctionCotangent<E>>, typename E::template deriv_type<1>, Order - 1>::type;

            FunctionCosecant(const E &fn) : fn(fn) { }

            template <size_t Order = 1>
            typename CosecantDerivativeType<E, Order>::type derivative(const Variable &var) const
            {
                if constexpr (Order == 0) return *this;
                else return (-csc(fn) * cot(fn) *
                    fn.template derivative<1>(var)).template derivative<Order-1>(var);
            }

            FunctionCosecant<E>& substituteInPlace(const std::vector<SubstituteEntry> &subs)
            {
                fn.substituteInPlace(subs);

                if (fn.isConcrete())
                {
                    _isConcrete = true;
                    _value = csc(fn.getValue());
                }

                return *this;
            }

            FunctionCosecant<E> substitute(const std::vector<SubstituteEntry> &subs) const
            {
                return FunctionCosecant<E>(*this).substituteInPlace(subs);
            }

            bool isConcrete() const { return _isConcrete; }

            double getValue() const { return _value; }

            friend std::ostream& operator<<(std::ostream &os, const FunctionCosecant<E> &cosecant)
            {
                if (cosecant._isConcrete) os << cosecant._value;
                else os << "csc(" << cosecant.fn << ")";

                return os;
            }

        private:
            E fn;

            bool _isConcrete = false;
            double _value;
    };

    template <typename E>
    FunctionCosecant<E> csc(FunctionExpression<E> const& fn)
    {
        return FunctionCosecant<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZTRIG_HH_

// vim: set ft=cpp.doxygen:
