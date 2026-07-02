# the submittable driver.
add_executable(p0004_main main.cpp)

# correctness: statement examples + edges + differential against the greedy oracle.
add_executable(p0004_tests tests.cpp)
add_test(NAME p0004 COMMAND p0004_tests)

# full solve at n = k = 2e5, measured end to end.
add_executable(p0004_bench bench.cpp)
