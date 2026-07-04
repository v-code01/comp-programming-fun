# the submittable driver.
add_executable(p0026_main main.cpp)

# correctness: classic hand shapes + enumerated edges + randomized differential
# against the floyd-warshall oracle on margin-clear graphs.
add_executable(p0026_tests tests.cpp)
add_test(NAME p0026 COMMAND p0026_tests)

# wall time of the O(n*m) detect at the worst legal bounds (n=300, m=90000).
add_executable(p0026_bench bench.cpp)
