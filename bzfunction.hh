#ifndef _BZFUNCTION_HH_
#define _BZFUNCTION_HH_

#include "bzvariable.hh"

#include <tuple>
#include <string>
#include <iostream>

namespace benzaiten
{
    struct FunctionImpl
    {
        virtual FunctionImpl* create() const = 0;
        virtual FunctionImpl* copy() const = 0;
        virtual bool isFunctionOfVar(const Variable &var) const = 0;
        virtual void print(std::ostream &os) const = 0;
        virtual void incrementDerivative(const Variable &var, size_t order = 1) = 0;
    };

    template <typename... Args>
    struct NullFunctionImpl : public FunctionImpl
    {
        NullFunctionImpl* create() const
        {
            return new NullFunctionImpl<Args...>;
        }

        NullFunctionImpl* copy() const
        {
            return new NullFunctionImpl<Args...>;
        }

        bool isFunctionOfVar(const Variable &var) const
        {
            // never a function of anything
            return false;
        }

        void incrementDerivative(const Variable &var, size_t order = 1)
        {
            // do nothing
        }

        void print(std::ostream &os) const
        {
            os << "<null>";
        }
    };

    template <typename... Args>
    struct AbstractFunctionImpl : public FunctionImpl
    {
        AbstractFunctionImpl* create() const
        {
            return new AbstractFunctionImpl<Args...>(name, args);
        }

        AbstractFunctionImpl* copy() const
        {
            return new AbstractFunctionImpl<Args...>(*this);
        }

        AbstractFunctionImpl(const std::string &name,
            const std::tuple<Args...> &args) :
                name(name), args(args)
        {
            for (size_t i = 0; i < sizeof...(Args); ++i) d[i] = 0;
        }

        bool isFunctionOfVar(const Variable &var) const
        {
            return isFunctionOfVarRecursive(var);
        }

        void incrementDerivative(const Variable &var, size_t order = 1)
        {
            incrementDerivativeRecursive(var, order);
        }

        void print(std::ostream &os) const
        {
            size_t order = totalDerivativeOrder();

            if (order > 0)
            {
                if (order > 1) os << "d^" << order;
                else os << "d";

                os << "(";
                os << name;
                os << "(";
                printArguments(os);
                os << ")";
                os << ")/";
                printDerivatives(os);
            }
            else
            {
                os << name;
                os << "(";
                printArguments(os);
                os << ")";
            }
        }

        private:
            /// Store the arguments to this function
            std::tuple<Args...> args;

            /// Unique name of this function
            std::string name;

            /// Derivatives and orders
            std::array<size_t, sizeof...(Args)> d;

            template <size_t I = 0>
            bool isFunctionOfVarRecursive(const Variable &var) const
            {
                if constexpr (I == sizeof...(Args))
                {
                    return false;
                }
                else if (std::get<I>(args).getName() == var.getName())
                {
                    return true;
                }
                else
                {
                    return isFunctionOfVarRecursive<I+1>(var);
                }
            }

            template <size_t I = 0>
            void incrementDerivativeRecursive(const Variable &var, size_t order = 1)
            {
                if constexpr (I == sizeof...(Args))
                {
                    return;
                }
                else if (std::get<I>(args).getName() == var.getName())
                {
                    d[I] += order;
                }
                else
                {
                    incrementDerivativeRecursive<I+1>(var, order);
                }
            }

            size_t totalDerivativeOrder() const
            {
                size_t order = 0;

                for (size_t i = 0; i < sizeof...(Args); ++i)
                {
                    order += d[i];
                }

                return order;
            }

            template <size_t I = 0>
            void printDerivatives(std::ostream &os) const
            {
                if constexpr (I == sizeof...(Args))
                {
                    return;
                }
                else
                {
                    if (d[I] > 0)
                    {
                        os << "d(" << std::get<I>(args).getName() << ")";
                        if (d[I] > 1) os << "^" << d[I];
                        os << " ";
                    }

                    printDerivatives<I+1>(os);
                }
            }

            template <size_t I = 0>
            void printArguments(std::ostream &os) const
            {
                if constexpr (I == sizeof...(Args))
                {
                    return;
                }
                else
                {
                    os << std::get<I>(args).getName();
                    if constexpr (I < sizeof...(Args) - 1) os << ", ";
                    printArguments<I+1>(os);
                }
            }
    };

    template <typename E>
    struct FunctionExpression
    {
    };

    struct Constant : public FunctionExpression<Constant>
    {
        template <size_t Order>
        using deriv_type = Constant;

        Constant(const double val) : value(val) { }

        template <size_t Order = 1>
        Constant& derivativeInPlace(const Variable &var)
        {
            value = 0;
            return *this;
        }

        template <size_t Order = 1>
        Constant derivative(const Variable &var) const
        {
            if constexpr (Order == 0) return *this;
            else return Constant(*this).derivativeInPlace<Order>(var);
        }

        friend std::ostream&
            operator<<(std::ostream &os, const Constant &cnst)
        {
            os << cnst.value;
            return os;
        }

        private:
            double value;
    };

    template <typename... Args>
    struct Function : public FunctionExpression<Function<Args...>>
    {
        template <size_t Order>
        using deriv_type = Function<Args...>;

        Function(const std::string &name, Args&... args)
        {
            impl = new AbstractFunctionImpl<Args...>(name, std::tuple<Args...>(args...));
        }

        Function(const Function<Args...> &other)
        {
            impl = other.impl->copy();
        }

        Function(Function<Args...> &&other)
        {
            impl = other.impl;
            other.impl = nullptr;
        }

        ~Function()
        {
            delete impl;
            impl = nullptr;
        }

        Function<Args...>& operator=(const Function<Args...> &other)
        {
            impl = other.impl->copy();
        }

        template <size_t Order = 1>
        Function<Args...>& derivativeInPlace(const Variable &var)
        {
            if (impl->isFunctionOfVar(var))
            {
                impl->incrementDerivative(var, Order);
            }
            else
            {
                delete impl;
                impl = new NullFunctionImpl<Args...>();
            }

            return *this;
        }

        template <size_t Order = 1>
        Function<Args...> derivative(const Variable &var) const
        {
            if constexpr (Order == 0) return *this;
            else return Function<Args...>(*this).derivativeInPlace<Order>(var);
        }

        friend std::ostream&
            operator<<(std::ostream &os, const Function<Args...> &fn)
        {
            fn.impl->print(os);
            return os;
        }

        private:
            /// Implementation pointer
            FunctionImpl *impl;
    };
}

#endif      // _BZFUNCTION_HH_

// vim: set ft=cpp.doxygen:
