# the submittable driver.
add_executable(p0036_main main.cpp)

# correctness: hand shapes + edges + differential vs the removal-reachability oracle.
add_executable(p0036_tests tests.cpp)
add_test(NAME p0036 COMMAND p0036_tests)

# n=2e5 m=3e5 worst case, measured.
add_executable(p0036_bench bench.cpp)
