# the submittable driver: stdin -> a_N mod p.
add_executable(p0037_main main.cpp)

# correctness: hand examples + edges + the differential vs both oracles.
add_executable(p0037_tests tests.cpp)
add_test(NAME p0037 COMMAND p0037_tests)

# BM + kitamasa wall time at m=2e4, order up to 1e4, N=1e18.
add_executable(p0037_bench bench.cpp)
