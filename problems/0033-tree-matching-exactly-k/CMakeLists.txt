# the submittable driver.
add_executable(p0033_main main.cpp)

# correctness: hand examples, edges, and aliens-vs-oracle over ALL K on thousands
# of small trees across shapes and weight regimes.
add_executable(p0033_tests tests.cpp)
add_test(NAME p0033 COMMAND p0033_tests)

# O(n log(n*maxW)) timing at n = 2e5, worst-case binary-search path.
add_executable(p0033_bench bench.cpp)
