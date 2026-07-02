# the submittable driver.
add_executable(p0008_main main.cpp)

# correctness: statement examples + edges + the differential vs the O(n^2 k) oracle.
add_executable(p0008_tests tests.cpp)
add_test(NAME p0008 COMMAND p0008_tests)

# wall time of the n=35000 k=50 worst case.
add_executable(p0008_bench bench.cpp)
