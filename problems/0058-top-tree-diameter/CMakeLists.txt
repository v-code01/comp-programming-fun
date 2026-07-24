# the submittable driver.
add_executable(p0058_main main.cpp)

# correctness: the machinery in isolation, then hand shapes + edges + differential
# vs the brute BFS/DFS diameter oracle.
add_executable(p0058_tests tests.cpp)
add_test(NAME p0058 COMMAND p0058_tests)

# n=1e5, q=2e5 adversarial mix, measured.
add_executable(p0058_bench bench.cpp)
