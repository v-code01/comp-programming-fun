# the submittable driver.
add_executable(p0016_main main.cpp)

# correctness: brute-verified examples + edges + the differential vs the subset-brute oracle.
add_executable(p0016_tests tests.cpp)
add_test(NAME p0016 COMMAND p0016_tests)

# wall time of the n=1e5, sum k=1e5 worst case.
add_executable(p0016_bench bench.cpp)
