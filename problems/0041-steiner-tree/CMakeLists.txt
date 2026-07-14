# the submittable driver.
add_executable(p0041_main main.cpp)

# correctness: hand shapes (forced hub, path, k=1, two terminals, all-terminals)
# + randomized differential against the 2^n subset + MST oracle, counting the
# cases where a Steiner point is what wins.
add_executable(p0041_tests tests.cpp)
add_test(NAME p0041 COMMAND p0041_tests)

# wall time of the subset DP at the worst legal bounds (n=100, k=10, complete).
add_executable(p0041_bench bench.cpp)
