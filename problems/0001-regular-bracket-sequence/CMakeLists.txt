# the submittable driver.
add_executable(p0001_main main.cpp)

# correctness: edges + the exhaustive grid against the search oracle.
add_executable(p0001_tests tests.cpp)
add_test(NAME p0001 COMMAND p0001_tests)

# throughput of the O(1) check.
add_executable(p0001_bench bench.cpp)
