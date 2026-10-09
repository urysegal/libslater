# Integrating libslater's four-center ERI into another service

Audience: a developer or coding agent adding STO electron repulsion integrals (ERIs) to a separate backend.
Read `README.md` first for dependencies, build commands and the C++ API. This document covers what you need
on top of that to run the library inside a long-lived service.

Everything below was checked on macOS (Apple Silicon, Homebrew) against the optimized build
(`build-release/`, `-O3`, no sanitizer). Items marked **unverified** were not checked.


## 1. What to take from this repo

You need the whole library, not individual files. The ERI code (`src/safouhi-4c.cpp` plus the f2c-translated
Fortran in `src/fortran-4c/`) depends on the rest of `src/`, `include/`, `precalculate/` and `libslater.h`.

Options, in order of preference:

1. **Git submodule / vendored copy + CMake `add_subdirectory`**. Configure with `-DSANITIZE=OFF
   -DCMAKE_BUILD_TYPE=Release`. Note that the top-level `CMakeLists.txt` also builds the example programs and the
   Catch2 test suite (fetched from GitHub at configure time), and it has no option to skip them.
2. **Build the shared library here and link against it.** You get `build-release/liblibslater.dylib` (`.so` on Linux).
   Link it with `-llibslater`, and also link `libf2c`, `spdlog` and `fmt` (see `README.md`). Set an rpath, or
   install the library.
3. **Run the `eri` executable as a subprocess** (`examples/eri.cpp`). This is the simplest option, and it isolates
   all the hazards in section 3. The cost is about 2 s of start-up per process (section 4). Instead of a one-shot
   CLI, consider a long-lived worker process that reads requests on stdin.

Build requirements that are easy to miss: Eigen **3.x** (Eigen 5 is rejected by CMake), and `libf2c`. On
Homebrew you also need `-DEigen3_DIR=/opt/homebrew/opt/eigen@3/share/eigen3/cmake` and
`-DCMAKE_EXE_LINKER_FLAGS=-L/opt/homebrew/lib -DCMAKE_SHARED_LINKER_FLAGS=-L/opt/homebrew/lib`.


## 2. The call

```cpp
#include "libslater.h"
using namespace slater;

// Once per process (or per worker):
STO_Integration_Engine factory;
std::map<integration_types, std::string> engines;      // empty = defaults
STO_Integrations *engine = factory.create(engines);
STO_Integration_Options options;
options.set(Quadrature_4C_algorithm_Parameter_Name, std::string("sdbar"));  // or "WGREP"
engine->init(options);

// Per integral:
STO_Basis_Function f1(STO_Basis_Function_Info(zeta1, {n1, l1, m1}), {Ax, Ay, Az});
// ... f2 @ B, f3 @ C, f4 @ D
double value = engine->electron_repulsion({f1, f2, f3, f4}).real();

// At shutdown:
delete engine;
```

### Semantics

- **Integral computed:** ⟨X₁(r₁−A) X₃(r₂−C) | 1/r₁₂ | X₂(r₁−B) X₄(r₂−D)⟩. Functions 1 and 2 belong to electron 1,
  and functions 3 and 4 to electron 2. In chemists' notation this is **(12|34)**. Map your backend's index
  convention to this carefully.
- **Units:** atomic units (bohr in, hartree out).
- **Normalization:** done internally. Pass only n, l, m, ζ.
- **Result:** real part of the returned `std::complex<double>`.
- **Real vs complex spherical harmonics:** **unverified**. If your basis uses real harmonics, establish the
  convention for m ≠ 0 against an independent code before relying on it.


## 3. Hazards for a service (important)

### 3.1 Not thread-safe
The Fortran-derived code keeps its working state in `static` variables throughout `src/fortran-4c/` (dozens
of files; `main.c` alone has ~70). Two concurrent `electron_repulsion` calls in the same process will corrupt
each other's results without any error. Either:
- serialize all calls with one global mutex (one integral at a time per process), or
- get parallelism from **multiple worker processes**, not threads.

### 3.2 Can terminate the process
f2c's `s_stop()` calls `exit()`. It is reachable from:
- `src/fortran-4c/gaunt.c`: when an internal Gaunt coefficient array is too small (high angular momentum), and
- `src/fortran-4c/zero_bsj.c`: when finding Bessel function roots fails during table initialization.

If this happens inside your backend process, the backend exits. Exceptions can't catch it. This is the main
reason to run the integrals in a separate worker process. The exact n/l limits that trigger it are
**unverified**. Test the largest quantum numbers your basis sets use before going to production.

### 3.3 Writes to stdout
Every call prints lines like `SDbar32_lrzv: ...`, `GREP-W32: ...` and a separator to **stdout** (Fortran
`WRITE` statements in `src/fortran-4c/main.c`). In a service, this will mix into whatever stdout is used for
(logs, or a stdin/stdout protocol with a worker). Either redirect the worker's stdout and use a different channel for
results, or remove the writes from `main.c`. If you remove them, keep that change in this repo so it isn't lost.

### 3.4 No input validation in optimized builds
`Quantum_Numbers::validate()` (0 ≤ l < n, |m| ≤ l) uses `assert`, which compiles away with `-DNDEBUG`.
Invalid input gives garbage or undefined behavior instead of an error. **Validate n, l, m and ζ > 0 in your
service before calling.**

### 3.5 Coincident centers
If all four functions are on the **same** center (a one-center integral), the result is ≈ 0, which is wrong.
The reference data itself has `3.36e-322` for these cases. Handle one-center integrals with a different method,
or reject them. Two-center cases where pairs coincide, like (aa|bb) with A ≠ B, **do** work; see section 5.


## 4. Performance (measured, optimized build, Apple Silicon)

| Situation | Time |
|---|---|
| First call in a process (builds quadrature and Bessel-root tables, ~24 MB static array) | ~2.2 s |
| Later calls, (1s 1s \| 2p₀ 2p₀) two-center | ~0.24 s each |
| Four distinct centers, 3d/2p/3p/2s mix | ~8.9 s (single process run, including the ~2 s start-up) |

Cost grows quickly with angular momentum. A full ERI tensor for a molecule is O(N⁴) of these calls, so plan
for batching, a work queue across worker processes, and caching (exploit the 8-fold permutational symmetry of
real integrals).


## 5. Known-good values for testing your integration

Use these as regression tests in the backend. Arguments are in `eri` CLI order: `n l m ζ x y z` × 4.

| Case | Inputs | Expected |
|---|---|---|
| 1s ×4, distinct centers on a line | `1 0 0 0.43 1.1 0 0` `1 0 0 1 2.2 0 0` `1 0 0 1 3.3 0 0` `1 0 0 1 4.4 0 0` | 0.2290246165 |
| with a 2p₋₁ | `1 0 0 0.43 1.1 0 0` `2 0 0 0.76 2.2 0 0` `2 1 -1 1.09 3.3 0 0` `1 0 0 1 4.4 0 0` | 0.07175464482 |
| (1s 1s \| 1s 1s), ζ=1, R=1.4 | `1 0 0 1 0 0 0` ×2, `1 0 0 1 0 0 1.4` ×2 | 0.503520932944 (= closed form 1/R − e^(−2R)(1/R + 11/8 + 3R/4 + R²/6)) |
| (1s 1s \| 2p₀ 2p₀), ζ=1.0/1.2, R=1.4 | `1 0 0 1.0 0 0 0` ×2, `2 1 0 1.2 0 0 1.4` ×2 | 0.474029188792 (no independent check) |

The first two come from `test/four-c-tests-data.cc`, which has about 340,000 more reference cases
(SDbar and GREP-W results) if you need a larger test set.


## 6. Suggested shape of the service

1. A small C++ worker binary (start from `examples/eri.cpp`) that initializes the engine once, then loops:
   read a request (four functions) → validate (3.4) → reject one-center cases (3.5) → compute → write the result.
   Put results on a channel separate from stdout, or remove the Fortran writes (3.3).
2. A pool of these workers managed by the backend, one integral at a time per worker (3.1). Restart workers that
   die (3.2), and treat a death as a failed request, not a backend crash.
3. Warm each worker at start-up with one dummy integral, to pay the ~2 s table build up front (section 4).
