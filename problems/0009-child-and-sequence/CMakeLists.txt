# the submittable driver: stdin -> stdout, fast IO, int64 sums.
add_executable(p0009_main main.cpp)

# correctness: both statement examples + enumerated edges + a differential
# sweep of thousands of random programs against the naive-array oracle.
add_executable(p0009_tests tests.cpp)
add_test(NAME p0009 COMMAND p0009_tests)

# amortized worst case at n = m = 1e5: adversarial set/mod/sum mix + mod storm.
add_executable(p0009_bench bench.cpp)
