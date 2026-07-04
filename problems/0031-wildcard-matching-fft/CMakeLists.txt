# the submittable driver.
add_executable(p0031_main main.cpp)

# correctness: hand examples, edges, the exhaustive small grid, the randomized
# differential, and the large-magnitude precision cases -- all vs the O(n*m) oracle.
add_executable(p0031_tests tests.cpp)
add_test(NAME p0031 COMMAND p0031_tests)

# worst-case wall clock of the two-FFT matcher at |T|=2e5, |P|=1e5.
add_executable(p0031_bench bench.cpp)
