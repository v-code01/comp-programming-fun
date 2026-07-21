# the submittable driver: stdin -> stdout, fast IO for n up to 1e6.
add_executable(p0050_main main.cpp)

# correctness: hand examples + edges + differential vs the O(n*V) DP oracle.
add_executable(p0050_tests tests.cpp)
add_test(NAME p0050 COMMAND p0050_tests)

# throughput of the O(n log n) slope-trick sweep at n=1e6.
add_executable(p0050_bench bench.cpp)
