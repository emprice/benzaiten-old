benzaiten: Symbolic differentiation in C++17
============================================

# Project overview

The use case for this project is as follows. Suppose you are solving a PDE numerically. You have an equation
with some known values (independent variables or functions of them) and unknown functions (what you are
solving for). Typically, when solving this kind of problem, you will have access to the values of
spatial and temporal derivatives of your unknown function, but, before you can use them, you need to expand
out all the terms in your equation.

You could do this with *Mathematica*, for example, but then you will need to transcribe those results into
code, which adds an element of human error to the process. You could use code generation software, like
SymPy; this removes the human element, but it can be very slow for large expressions. It also adds an extra
step to compiling your program: Generate code in Python, then compile, then run. This can cause problems of
its own if the dependency chain is not set up properly.

What I would like to do in this situation is write my equations down in pure C++ and have all those functions
and derivatives expand themselves, all in one program. This is the purpose of benzaiten. Benzaiten does not 
take numerical derivatives. Rather, it takes analytic derivatives, propagating derivatives of unknown 
functions, and substitutes values at the end. Lots of work is done directly in the compiling step, because
benzaiten is a pure C++17 template library. Expressions should be reusable, perfect for working in loops.

# How do I use it?

Check out [bztest.cc](https://github.com/emprice/benzaiten/blob/master/bztest.cc) for an example of using
benzaiten to evaluate simple and complex expressions in pure C++17. I chose C++17 to take advantage
of some newer features; it is the only requirement for compiling code against this library, as there
are no external dependencies.

# Why the name?

[Benzaiten](https://en.wikipedia.org/wiki/Benzaiten) is the Japanese goddess of "everything that flows",
including water. This seemed an appropriate name for software that was developed to be part of a
larger hydrodynamics code.
