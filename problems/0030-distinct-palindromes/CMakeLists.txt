# the submittable driver.
add_executable(p0030_main main.cpp)

# correctness: hand examples + edges + exhaustive grids + randomized differential
# against the O(n^2) set oracle.
add_executable(p0030_tests tests.cpp)
add_test(NAME p0030 COMMAND p0030_tests)

# wall-clock of the eertree build at |s|=1e6.
add_executable(p0030_bench bench.cpp)
