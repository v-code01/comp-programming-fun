# the submittable driver.
add_executable(p0038_main main.cpp)

# correctness: hand-checked shapes + randomized differential against the
# independent per-pair Edmonds-Karp oracle, including disconnected 0-cuts.
add_executable(p0038_tests tests.cpp)
add_test(NAME p0038 COMMAND p0038_tests)

# wall time to build the tree at the worst bounds (n=150, complete graph).
add_executable(p0038_bench bench.cpp)
