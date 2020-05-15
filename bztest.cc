#include "benzaiten.hh"

using namespace benzaiten;

int main(int argc, char **argv)
{
    Variable t("t", Temporal), x("x", Spatial), y("y", Spatial);
    Function f("f", t, x, y);
    Function g("g", x, t);
    Function h("h", t);

    // testing derivative and printing
    std::cout << "<<< testing multi derivative >>>" << std::endl;
    auto df = f.derivative(x).derivative<2>(y);
    std::cout << df << std::endl << std::endl;

    // testing addition
    std::cout << "<<< testing addition >>>" << std::endl;
    auto sum1 = f + g;
    std::cout << sum1.derivative(x).derivative(y) << std::endl;
    auto sum2 = f + 1;
    std::cout << sum2.derivative(x).derivative(y) << std::endl;
    auto sum3 = 1 + f;
    std::cout << sum3.derivative(x).derivative(y) << std::endl << std::endl;

    // testing subtraction
    std::cout << "<<< testing subtraction >>>" << std::endl;
    auto diff1 = f - g;
    std::cout << diff1.derivative(x).derivative(y) << std::endl;
    auto diff2 = f - 1;
    std::cout << diff2.derivative(x).derivative(y) << std::endl;
    auto diff3 = 1 - f;
    std::cout << diff3.derivative(x).derivative(y) << std::endl << std::endl;

    // testing multiplication
    std::cout << "<<< testing multiplication >>>" << std::endl;
    auto prod1 = f * g;
    std::cout << prod1.derivative<2>(x) << std::endl;
    auto prod2 = f * 2;
    std::cout << prod2.derivative<2>(x) << std::endl;
    auto prod3 = 2 * f;
    std::cout << prod3.derivative<2>(x) << std::endl << std::endl;

    // testing division
    std::cout << "<<< testing division >>>" << std::endl;
    auto quo1 = f / g;
    std::cout << quo1.derivative(x) << std::endl;
    auto quo2 = f / 3.;
    std::cout << quo2.derivative(x) << std::endl;
    auto quo3 = 3. / f;
    std::cout << quo3.derivative(x) << std::endl << std::endl;

    // testing powers
    std::cout << "<<< testing function power >>>" << std::endl;
    auto pwr = f ^ g;
    std::cout << pwr.derivative(x) << std::endl << std::endl;

    std::cout << "<<< testing square root >>>" << std::endl;
    auto rt = sqrt(f);
    std::cout << rt.derivative<2>(x) << std::endl << std::endl;

    // testing log and exp
    std::cout << "<<< testing log and exp >>>" << std::endl;
    auto lg = log(f);
    std::cout << lg.derivative(x) << std::endl;

    auto ex = exp(f);
    std::cout << ex.derivative(x) << std::endl << std::endl;

    // testing trigonometry
    std::cout << "<<< testing trig >>>" << std::endl;

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
    std::cout << cosecant.derivative(x) << std::endl << std::endl;;

    // testing hyperbolic trigonometry
    std::cout << "<<< testing hyperbolic trig >>>" << std::endl;

    auto sinh_ = sinh(f);
    std::cout << sinh_.derivative<2>(x) << std::endl;

    auto cosh_ = cosh(f);
    std::cout << cosh_.derivative<2>(x) << std::endl;

    auto tanh_ = tanh(f);
    std::cout << tanh_.derivative(x) << std::endl;

    auto coth_ = coth(f);
    std::cout << coth_.derivative(x) << std::endl;

    auto sech_ = sech(f);
    std::cout << sech_.derivative(x) << std::endl;

    auto csch_ = csch(f);
    std::cout << csch_.derivative(x) << std::endl << std::endl;

    // testing variable expressions
    std::cout << "<<< testing variable expressions >>>" << std::endl;
    auto expr = (f * (x ^ 2));
    std::cout << expr.derivative(x) << std::endl;

    auto expr2 = h * log(x);
    std::cout << expr2.derivative(x) << std::endl;

    auto expr3 = (3. / sqrt(1. / (x ^ 3))) * f * sqrt(x);
    std::cout << expr3.derivative(x) << std::endl;

    auto expr4 = x * (t ^ 2);
    std::cout << expr4.derivative(x).derivative(x) << std::endl << std::endl;

    // testing substitute
    std::cout << "<<< testing substitution >>>" << std::endl;
    std::vector<SubstituteEntry> subs;
    subs.push_back(SubstituteEntry("f", 1, { }));
    subs.push_back(SubstituteEntry("g", 2, { }));
    subs.push_back(SubstituteEntry("f", 3, { { "x", 1 } }));
    subs.push_back(SubstituteEntry("g", 4, { { "x", 1 } }));
    subs.push_back(SubstituteEntry("x", 5, { }));
    subs.push_back(SubstituteEntry("t", 6, { }));

    auto test1 = (f * csc(g)).derivative(x);
    std::cout << test1 << " = " << test1.substitute(subs) << std::endl;

    auto test2 = (f / csch(g)).derivative(x);
    std::cout << test2 << " = " << test2.substitute(subs) << std::endl;

    auto test3 = (f * (x ^ 2)).derivative(x);
    std::cout << test3 << " = " << test3.substitute(subs) << std::endl;

    auto test4 = ((3. / sqrt(1. / (x ^ 3))) * f * sqrt(x)).derivative(x);
    std::cout << test4 << " = " << test4.substitute(subs) << std::endl;

    auto test5 = (x * (x * (t ^ 2)).derivative(x)).derivative(x);
    std::cout << test5 << " = " << test5.substitute(subs) << std::endl << std::endl;

    auto test6 = (1. / x) * (x * f).derivative(x);
    std::cout << test6 << " = " << test6.substitute(subs) << std::endl << std::endl;

    return 0;
}

// vim: set ft=cpp.doxygen:
