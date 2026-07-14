# the submittable driver.
add_executable(p0042_main main.cpp)

# correctness: the ballot pin + edges + the compatibility hypothesis + a
# differential against the brute path enumeration.
add_executable(p0042_tests tests.cpp)
add_test(NAME p0042 COMMAND p0042_tests)

# k=100, coords 2e5, dense M -- the worst case the constraints allow.
add_executable(p0042_bench bench.cpp)
