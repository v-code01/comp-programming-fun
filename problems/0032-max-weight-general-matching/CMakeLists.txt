# the submittable driver.
add_executable(p0032_main main.cpp)

# correctness: hand shapes + differential vs the 2^n subset-DP oracle.
add_executable(p0032_tests tests.cpp)
add_test(NAME p0032 COMMAND p0032_tests)

# n=500 complete-graph worst case, measured.
add_executable(p0032_bench bench.cpp)
