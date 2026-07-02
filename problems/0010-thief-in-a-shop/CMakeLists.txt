# the submittable driver.
add_executable(p0010_main main.cpp)

# correctness: statement examples, edges, and the differential vs the set oracle.
add_executable(p0010_tests tests.cpp)
add_test(NAME p0010 COMMAND p0010_tests)

# worst-case wall clock of the O(V*D) min-coin dp.
add_executable(p0010_bench bench.cpp)
