# the submittable driver.
add_executable(p0003_main main.cpp)

# correctness: statement examples + edges + differential vs the brute oracle.
add_executable(p0003_tests tests.cpp)
add_test(NAME p0003 COMMAND p0003_tests)

# full-solve wall time at the n=q=5000 worst case.
add_executable(p0003_bench bench.cpp)
