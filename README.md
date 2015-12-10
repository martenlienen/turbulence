# Turbulent Flow Simulation

## Libraries

We require the following libraries

- **openmpi** version *1.10.1*
- **petsc** version *3.5.4*
- **gtest** version *1.6.0*

## Setup

For development you should install the following tools

* `clang`: We use `clang`s `clang-format` to keep our coding style consistent.
* `cmake`: We use `cmake` to generate a nicer `Makefile` than the one we got
  from the tutors.
* `gtest`: GTest is used for unit testing. Run the tests with `./tests`.

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

## Connecting to the cluster

```sh
# Read the usage information of the cluster connection script
scripts/cluster/connect -h

# Connect to cluster
scripts/cluster/connect ...

# On the cluster

# Clone the repository
git clone https://github.com/cqql/turbulence

# Move into it
cd turbulence

# Load necessary modules (PETSc etc.)
#
# *Note the . at the beginning of the line!* This makes the script modify the
# environment of your current shell instead of creating a subshell for the
# execution of the script.
. scripts/cluster/load-modules

# Compile the simulation
cmake -Dtests=OFF .
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
