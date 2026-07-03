# the submittable driver.
add_executable(p0020_main main.cpp)

# correctness: statement examples, edges, and the differential vs the brute oracle.
add_executable(p0020_tests tests.cpp)
add_test(NAME p0020 COMMAND p0020_tests)

# dense worst-case wall clock of the min-cost max-flow.
add_executable(p0020_bench bench.cpp)
