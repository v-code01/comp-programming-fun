# the submittable driver.
add_executable(p0055_main main.cpp)

# correctness: hand examples + edges, differential vs trial division (n<=1e12),
# and the multiply-back + primality self-check on the full 1e18 range.
add_executable(p0055_tests tests.cpp)
add_test(NAME p0055 COMMAND p0055_tests)

# wall time to factor 500 hard ~1e18 semiprimes, plus steady-state throughput.
add_executable(p0055_bench bench.cpp)
