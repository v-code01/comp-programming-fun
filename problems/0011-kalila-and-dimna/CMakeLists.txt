# the submittable driver.
add_executable(p0011_main main.cpp)

# correctness: statement examples + edges + differential against the O(n^2) oracle.
add_executable(p0011_tests tests.cpp)
add_test(NAME p0011 COMMAND p0011_tests)

# full solve at n = 1e5, measured end to end.
add_executable(p0011_bench bench.cpp)
