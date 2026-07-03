# the submittable driver: stdin -> stdout, fast IO, iterative Euler tour.
add_executable(p0017_main main.cpp)

# correctness: the brute-verified 620E sample + enumerated edges (single node,
# path, star, whole-tree repaint, nested repaints, all-same-color, color ceiling)
# + a differential sweep of thousands of random programs on random trees against
# the naive-array oracle.
add_executable(p0017_tests tests.cpp)
add_test(NAME p0017 COMMAND p0017_tests)

# worst case at n = m = 4e5: deep path tree + bushy random tree.
add_executable(p0017_bench bench.cpp)
