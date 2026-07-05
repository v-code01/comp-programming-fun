# the submittable driver.
add_executable(p0034_main main.cpp)

# correctness: OEIS cross-checks + edges + differential vs the brute-DFS oracle.
add_executable(p0034_tests tests.cpp)
add_test(NAME p0034 COMMAND p0034_tests)

# full 12x12 and 10x12 worst cases, measured.
add_executable(p0034_bench bench.cpp)
