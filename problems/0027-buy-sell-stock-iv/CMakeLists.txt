# the submittable driver.
add_executable(p0027_main main.cpp)

# correctness: hand examples, edges, oracle-vs-brute, aliens-vs-oracle over all K.
add_executable(p0027_tests tests.cpp)
add_test(NAME p0027 COMMAND p0027_tests)

# O(n log maxPrice) timing at n = 1e6, worst-case binary-search path.
add_executable(p0027_bench bench.cpp)
