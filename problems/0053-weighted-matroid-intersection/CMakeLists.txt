# the submittable driver.
add_executable(p0053_main main.cpp)

# correctness: hand shapes + edges + differential vs the 2^m brute oracle.
add_executable(p0053_tests tests.cpp)
add_test(NAME p0053 COMMAND p0053_tests)

# n=m=200 worst case, measured.
add_executable(p0053_bench bench.cpp)
