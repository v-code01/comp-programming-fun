# the submittable driver.
add_executable(p0021_main main.cpp)

# correctness: hand-checked edges + randomized differential against the 2^n
# subset-closure oracle, including prerequisite cycles and negative profits.
add_executable(p0021_tests tests.cpp)
add_test(NAME p0021 COMMAND p0021_tests)

# wall time of the max-flow solve at the worst legal bounds (n=5000, m=20000).
add_executable(p0021_bench bench.cpp)
