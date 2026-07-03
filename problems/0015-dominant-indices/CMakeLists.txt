# the submittable driver -- fast IO, iterative long-path decomposition over stdin.
add_executable(p0015_main main.cpp)

# correctness: statement shapes + edges + structured sweep + randomized
# differential against the O(n^2) histogram oracle.
add_executable(p0015_tests tests.cpp)
add_test(NAME p0015 COMMAND p0015_tests)

# solve wall time at n=1e6, including the pathological path shape.
add_executable(p0015_bench bench.cpp)
