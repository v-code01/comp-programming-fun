# the submittable driver -- reads the tree over stdin, dsu-on-tree, answers out.
add_executable(p0013_main main.cpp)

# correctness: statement example + enumerated edges (single/path/star/all-same/
# all-distinct/tie) + randomized differential against the O(n^2) subtree oracle.
add_executable(p0013_tests tests.cpp)
add_test(NAME p0013 COMMAND p0013_tests)

# runtime across worst-case shapes at n=1e5, including a depth-1e5 path.
add_executable(p0013_bench bench.cpp)
