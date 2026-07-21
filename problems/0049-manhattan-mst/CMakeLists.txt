# the submittable driver.
add_executable(p0049_main main.cpp)

# correctness: hand cases + edges + the differential against the O(n^2) oracle.
add_executable(p0049_tests tests.cpp)
add_test(NAME p0049 COMMAND p0049_tests)

# throughput of the octant sweep at n=1e5.
add_executable(p0049_bench bench.cpp)
