# comp-programming-fun

i love solving these in my free time :)

every problem in here gets the full treatment -- not the first solution that
passes, the best one that exists. optimal complexity, and a note on why it can't
go faster. checked against a brute-force oracle on thousands of random cases.
benchmarked, so the numbers are real.

one folder per problem. the solution, the naive version i diff it against, the
tests, the benchmark.

that's it. have fun.

## build

```sh
cmake -B build -G Ninja
cmake --build build
```

tests run under sanitizers -- use the `Sanitize` build type for the leak and ub
checks:

```sh
cmake -B build-san -G Ninja -DCMAKE_BUILD_TYPE=Sanitize
cmake --build build-san
```

## layout

```
problems/NNNN-slug/
  solution.hpp    the fast one
  reference.hpp   the dumb one -- obviously correct, what i diff against
  main.cpp        stdin -> solution -> stdout
  tests.cpp       edge cases + differential (or exhaustive when the space is small)
  bench.cpp       measured ns/op
```

## other languages

most of these are c++. a few are in go or python where that was the better tool
-- same deal, same rigor: a fast solution, a dumb reference, a differential
harness, a benchmark. the c++ build ignores those folders (no CMakeLists inside),
so they never touch the cmake step.

```sh
go test cpf/...                          # every go problem, via the go.work
python3 problems/<slug>/test_solution.py # a python problem's differential
```
