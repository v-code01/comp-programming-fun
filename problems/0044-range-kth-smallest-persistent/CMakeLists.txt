# the submittable driver.
add_executable(p0044_main main.cpp)

# correctness: hand examples + edges + randomized differential vs the sort oracle.
add_executable(p0044_tests tests.cpp)
add_test(NAME p0044 COMMAND p0044_tests)

# build + query throughput at n = q = 1e5.
add_executable(p0044_bench bench.cpp)
