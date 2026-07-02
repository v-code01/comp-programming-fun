# the submittable driver.
add_executable(p0007_main main.cpp)

# correctness: statement examples + edges + the differential vs the 2^n brute.
add_executable(p0007_tests tests.cpp)
add_test(NAME p0007 COMMAND p0007_tests)

# throughput of the SOS solve at n = 1e6, 20-bit values.
add_executable(p0007_bench bench.cpp)
