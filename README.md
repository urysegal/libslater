# libslater
C++ library to calculate various ab initio chemistry integrals using Slater basis sets

This code implements in C++ the following works:

[1] Slevinsky, R.M., Safouhi, H. Compact formulae for three-center nuclear attraction integrals over exponential type functions. J Math Chem 60, 1337–1355 (2022).


[2] H.H.H. Homeier, E.J. Weniger, E.O. Steinborn,
Programs for the evaluation of overlap integrals with B functions, 
Computer Physics Communications, Volume 72, Issues 2–3, 1992, Pages 269-287, ISSN 0010-4655


[3] On the evaluation of overlap integrals with exponential‐type basis functions, 
HHH Homeier, EO Steinborn,  International journal of quantum chemistry 42 (4), 761-778

The four-center electron repulsion integrals are computed by an f2c translation of
H. Safouhi's Fortran code (`src/fortran-4c/`), using the S\bar{D} and GREP-W methods.


If you're embedding the ERI code in a service or another backend, also read
[`INTEGRATING.md`](INTEGRATING.md). It covers thread safety, process-exit hazards, stdout output and performance.


## Dependencies

| Package | Notes |
|---------|-------|
| CMake ≥ 3.16 | |
| C++20 compiler | Apple clang works |
| Boost ≥ 1.71 | header-only parts (math, quadrature) |
| Eigen **3.x** | Eigen 5 is rejected by `find_package(Eigen3 3.3)` |
| spdlog, fmt | logging |
| libf2c | runtime for the translated Fortran code |
| Catch2 v3 | fetched automatically by CMake for the tests |

On macOS with Homebrew:

```
brew install cmake boost eigen@3 spdlog fmt libf2c
```

Homebrew's default `eigen` formula is Eigen 5, which won't work. Install `eigen@3` and point CMake at it,
as shown below.


## Building

There are two build flavors, each in its own build directory.

### Debug (default): AddressSanitizer, `-g`, no optimization

```
cmake -S . -B build \
  -DEigen3_DIR=/opt/homebrew/opt/eigen@3/share/eigen3/cmake \
  -DCMAKE_EXE_LINKER_FLAGS=-L/opt/homebrew/lib \
  -DCMAKE_SHARED_LINKER_FLAGS=-L/opt/homebrew/lib
cmake --build build -j8
```

### Optimized: `-O3 -DNDEBUG`, no sanitizer

```
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF \
  -DEigen3_DIR=/opt/homebrew/opt/eigen@3/share/eigen3/cmake \
  -DCMAKE_EXE_LINKER_FLAGS=-L/opt/homebrew/lib \
  -DCMAKE_SHARED_LINKER_FLAGS=-L/opt/homebrew/lib
cmake --build build-release -j8
```

Use the optimized build for real calculations. The debug build is much slower because of AddressSanitizer.

Notes:
- `SANITIZE` (default `ON`) controls `-fsanitize=address -g`. Setting only `CMAKE_BUILD_TYPE=Release`
  is not enough, because the sanitizer flags are added per-target.
- The `-L/opt/homebrew/lib` flags are needed so the linker finds `libf2c`. On Linux with system packages
  you can usually leave out the Homebrew-specific flags.
- `-Wall -Werror` is on for the C++ code (it's relaxed for `src/fortran-4c/`), so new warnings break the build.
- Each build directory produces `liblibslater.dylib` (or `.so`), the example programs, and `test/test-all`.


## Running the tests

```
cd build-release && ctest -j8 --output-on-failure
```

You can also run `build-release/test/test-all` directly (it's a Catch2 binary, so you can pass a test name or `[tag]`).
The four-center test checks about 340,000 reference values from `test/four-c-tests-data.cc`. Run it in
the optimized build: under AddressSanitizer it takes a very long time.


## Computing a four-center ERI

### From the command line

The `eri` example program takes four basis functions, each given as `n l m zeta x y z`:

```
build-release/eri  n1 l1 m1 zeta1 Ax Ay Az   n2 l2 m2 zeta2 Bx By Bz \
                   n3 l3 m3 zeta3 Cx Cy Cz   n4 l4 m4 zeta4 Dx Dy Dz
```

Example:

```
$ build-release/eri  1 0 0 0.43 1.1 0 0   2 0 0 0.76 2.2 0 0   2 1 -1 1.09 3.3 0 0   1 0 0 1 4.4 0 0
...
0.0717546448169
```

The last line is the integral. The Fortran code also prints its own diagnostic lines
(`SDbar32_lrzv: ...`, `GREP-W32: ...`) to stdout on every call, so use `tail -1` if you're parsing the output.

### From C++

```cpp
#include "libslater.h"
using namespace slater;

STO_Basis_Function f1(STO_Basis_Function_Info(zeta1, {n1, l1, m1}), {Ax, Ay, Az});
STO_Basis_Function f2(STO_Basis_Function_Info(zeta2, {n2, l2, m2}), {Bx, By, Bz});
STO_Basis_Function f3(STO_Basis_Function_Info(zeta3, {n3, l3, m3}), {Cx, Cy, Cz});
STO_Basis_Function f4(STO_Basis_Function_Info(zeta4, {n4, l4, m4}), {Dx, Dy, Dz});

STO_Integration_Engine factory;
std::map<integration_types, std::string> engines;   // empty map = default implementations
STO_Integrations *engine = factory.create(engines);

STO_Integration_Options options;
options.set(Quadrature_4C_algorithm_Parameter_Name, std::string("sdbar"));  // or "WGREP"
engine->init(options);

double eri = engine->electron_repulsion({f1, f2, f3, f4}).real();

delete engine;
```

See `examples/eri.cpp` for a complete program. To compile your own program against the library:

```
c++ -std=c++20 -O2 -I<libslater> -I/opt/homebrew/include my_prog.cpp \
    -L<libslater>/build-release -llibslater -Wl,-rpath,<libslater>/build-release
```

### Conventions

- **Electron assignment.** The integral computed is

  ⟨ X₁(r₁−A) X₃(r₂−C) | 1/|r₁−r₂| | X₂(r₁−B) X₄(r₂−D) ⟩

  Functions 1 and 2 belong to electron 1, and functions 3 and 4 to electron 2. In chemists' notation that is (12|34).
- **Basis functions.** These are STOs with radial part N·r^(n−1)·e^(−ζr) times a spherical harmonic Y_l^m.
  The library normalizes each function itself, so you only pass n, l, m and ζ.
- **Quantum numbers.** They must satisfy 0 ≤ l < n and |m| ≤ l. This is only checked with `assert`, so the
  optimized build accepts invalid input silently.
- **Algorithm.** Both S\bar{D} (default) and GREP-W are computed on every call, and the option picks
  which one is returned. They normally agree to about 10 digits.
- **Return type.** The result is `std::complex<double>`. Use `.real()`.

### Limitations

- **Coincident centers.** The code doesn't handle the case where all four centers coincide (one-center
  integrals): it returns ≈ 0 (the reference data contains `3.36e-322` for these cases). Use a different
  method for one-center integrals.
- Not verified here: whether nonzero `m` means complex or real spherical harmonics. Check this before relying on
  results with m ≠ 0 in a real-harmonic basis.


## Other integrals

`STO_Integrations` also provides:

- `overlap({f, g})`
- `kinetic({f, g})`
- `nuclear_attraction({f, g}, nucleus_center)`

See `libslater.h` for the full API. The programs in `examples/` show usage of each integral type:

| Program | What it does |
|---------|--------------|
| `eri` | Four-center ERI with arbitrary n, l, m, ζ and centers (see above) |
| `simple4c_any` | Four-center ERI for 1s functions only, with ζ and centers given as `-a1 ζ -c1 x y z ...` (used by `4c_any.sh`) |
| `simple4c` | Four-center ERI with hardcoded inputs |
| `test4c <id> [to-id]` | Recomputes reference cases from the 4-center test data and compares |
| `test3c` | Same idea for three-center nuclear attraction |
| `simple`, `simple2`, `simple3` | Two-center overlap (and kinetic, in `simple`) examples |
