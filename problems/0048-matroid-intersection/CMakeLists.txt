# the submittable driver.
add_executable(p0048_main main.cpp)

# correctness: hand shapes + edges + differential vs the 2^m brute oracle.
add_executable(p0048_tests tests.cpp)
add_test(NAME p0048 COMMAND p0048_tests)

# n=m=300 worst case, measured.
add_executable(p0048_bench bench.cpp)
