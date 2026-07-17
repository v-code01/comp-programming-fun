# the submittable driver: stdin -> stdout, signed fast IO, int64 sums.
add_executable(p0045_main main.cpp)

# correctness: hand examples (a chmin that beats, add-then-chmin) + enumerated
# edges (no-op cap, flatten, single element, all-equal, negative adds) + a
# differential sweep of thousands of beats-heavy programs against the array oracle.
add_executable(p0045_tests tests.cpp)
add_test(NAME p0045 COMMAND p0045_tests)

# amortized worst case at n = q = 1e5: chmin storm + descending cap + mixed.
add_executable(p0045_bench bench.cpp)
