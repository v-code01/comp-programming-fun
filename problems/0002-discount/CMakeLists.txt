# the submittable driver.
add_executable(p0002_main main.cpp)

# correctness: statement example + edges + differential against the O(n^2) oracle.
add_executable(p0002_tests tests.cpp)
add_test(NAME p0002 COMMAND p0002_tests)

# throughput of the full sort-and-answer solve.
add_executable(p0002_bench bench.cpp)
