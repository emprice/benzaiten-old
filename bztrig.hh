#ifndef _BZTRIG_HH_
#define _BZTRIG_HH_

#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzneg.hh"

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

            friend std::ostream& operator<<(std::ostream &os, const FunctionSine<E> &sine)
            {
                os << "sin(" << sine.fn << ")";
                return os;
            }

        private:
            const E fn;
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

            friend std::ostream& operator<<(std::ostream &os, const FunctionCosine<E> &cosine)
            {
                os << "cos(" << cosine.fn << ")";
                return os;
            }

        private:
            const E fn;
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

            friend std::ostream& operator<<(std::ostream &os, const FunctionTangent<E> &tngt)
            {
                os << "tan(" << tngt.fn << ")";
                return os;
            }

        private:
            const E fn;
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

            friend std::ostream& operator<<(std::ostream &os, const FunctionCotangent<E> &cotan)
            {
                os << "cot(" << cotan.fn << ")";
                return os;
            }

        private:
            const E fn;
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

            friend std::ostream& operator<<(std::ostream &os, const FunctionSecant<E> &secant)
            {
                os << "sec(" << secant.fn << ")";
                return os;
            }

        private:
            const E fn;
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

            friend std::ostream& operator<<(std::ostream &os, const FunctionCosecant<E> &cosecant)
            {
                os << "csc(" << cosecant.fn << ")";
                return os;
            }

        private:
            const E fn;
    };

    template <typename E>
    FunctionCosecant<E> csc(FunctionExpression<E> const& fn)
    {
        return FunctionCosecant<E>(static_cast<E const&>(fn));
    }
}

#endif      // _BZTRIG_HH_

// vim: set ft=cpp.doxygen:
