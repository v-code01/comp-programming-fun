# the submittable driver.
add_executable(p0006_main main.cpp)

# correctness: statement examples + edges + exhaustive small vs deletion brute
# force + randomized differential vs the last-char oracle.
add_executable(p0006_tests tests.cpp)
add_test(NAME p0006 COMMAND p0006_tests)

# wall clock of the O(n^3) at the n=500 ceiling.
add_executable(p0006_bench bench.cpp)
