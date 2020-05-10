#include "bzvariable.hh"
#include "bzfunction.hh"
#include "bzsum.hh"
#include "bzdifference.hh"
#include "bzproduct.hh"
#include "bzpower.hh"
#include "bzlog.hh"
#include "bztrig.hh"

using namespace benzaiten;

int main(int argc, char **argv)
{
    Variable t("t", Temporal), x("x", Spatial), y("y", Spatial);
    Function f("f", t, x, y);
    Function g("g", x, t);
    Function h("h", x, t, y);

    // testing derivative and printing
    auto df = f.derivative(x).derivative<2>(y);
    std::cout << df << std::endl;

    // testing addition
    auto sum = f + g;
    std::cout << sum.derivative(x).derivative(y) << std::endl;

    // testing subtraction
    auto diff = f - g;
    std::cout << diff.derivative(x).derivative(y) << std::endl;

    // testing multiplication
    auto prod = f * g;
    std::cout << prod.derivative<2>(x) << std::endl;

    // testing division
    auto quo = f / g;
    std::cout << quo.derivative(x) << std::endl;

    // testing trigonometry
    auto sine = sin(f);
    std::cout << sine.derivative<2>(x) << std::endl;

    auto cosine = cos(f);
    std::cout << cosine.derivative<2>(x) << std::endl;

    auto tangent = tan(f);
    std::cout << tangent.derivative(x) << std::endl;

    auto cotangent = cot(f);
    std::cout << cotangent.derivative(x) << std::endl;

    auto secant = sec(f);
    std::cout << secant.derivative(x) << std::endl;

    auto cosecant = csc(f);
    std::cout << cosecant.derivative(x) << std::endl;

    // testing variable expressions
    auto expr = (f * (x ^ Constant(2)));
    std::cout << (x ^ Constant(2)).derivative(x) << std::endl;

    // testing substitute
    std::vector<SubstituteEntry> subs;
    subs.push_back(SubstituteEntry("f", 1, { }));
    subs.push_back(SubstituteEntry("g", 2, { }));
    subs.push_back(SubstituteEntry("f", 3, { { "x", 1 } }));
    subs.push_back(SubstituteEntry("g", 4, { { "x", 1 } }));

    auto result = (f * csc(g)).derivative(x).substitute(subs);
    std::cout << result << std::endl;

    return 0;
}

// vim: set ft=cpp.doxygen:
