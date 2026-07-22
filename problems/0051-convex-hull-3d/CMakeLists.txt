# the submittable driver.
add_executable(p0051_main main.cpp)

# correctness: hand cases + edges + the differential against the O(n^4) oracle.
add_executable(p0051_tests tests.cpp)
add_test(NAME p0051 COMMAND p0051_tests)

# throughput of the incremental hull at n=2000.
add_executable(p0051_bench bench.cpp)
