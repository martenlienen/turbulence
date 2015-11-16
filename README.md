# Turbulent Flow Simulation

## Libraries

We require the following libraries

- **openmpi** version *1.10.1*
- **petsc** version *3.5.4*

## Setup

For development you should install the following tools

* `clang`: We use `clang`s `clang-format` to keep our coding style consistent.
* `cmake`: We use `cmake` to generate a nicer `Makefile` than the one we got
  from the tutors.

```sh
# Run these once at the start of your session or put it into your
# .bashrc/.bash_profile. We do not want these values in our Makefile/
# CMakeLists.txt because they will be different for everyone and would thus
# generate lots of useless change sets.
export PETSC_DIR=<..>
export PETSC_ARCH=<..>

# Run this once in the beginning to generate a nice Makefile. You only need to
# run this again when you modified the CMakeLists.txt, for example when you add
# a new .cpp file.
cmake .

# Recompiling the program just needs a simple make. This make also has proper
# dependency tracking so that you should basically never need to use make clean.
make
```

## Style Guide

We use
[google's C++ style guide](http://google.github.io/styleguide/cppguide.html). It
is a reasonable guideline, used by a large organization to great success and at
the same time it is an interesting read and can maybe even teach you something
about C++.

You can use the script in `scripts/format` to automatically reformat all source
files. Alternatively try one of the integrations of `clang-format` into various
editors/IDEs.
