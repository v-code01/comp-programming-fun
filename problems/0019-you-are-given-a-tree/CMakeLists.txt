# the submittable driver -- fast IO, iterative greedy + sqrt-on-the-answer.
add_executable(p0019_main main.cpp)

# correctness: the statement example + edges + structured sweep + randomized
# differential against the independent exponential set-packing oracle across all k.
add_executable(p0019_tests tests.cpp)
add_test(NAME p0019 COMMAND p0019_tests)

# full n-answer solve wall time at n=1e5 across path/star/random shapes.
add_executable(p0019_bench bench.cpp)
