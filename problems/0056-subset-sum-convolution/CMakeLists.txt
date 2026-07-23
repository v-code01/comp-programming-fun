# the submittable driver.
add_executable(p0056_main main.cpp)

# correctness: hand examples + edges + the differential vs the O(3^k) submask oracle.
add_executable(p0056_tests tests.cpp)
add_test(NAME p0056 COMMAND p0056_tests)

# throughput of the ranked zeta/mobius solve at k=20.
add_executable(p0056_bench bench.cpp)
