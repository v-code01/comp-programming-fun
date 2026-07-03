# the submittable driver -- centroid tree over stdin, fast IO.
add_executable(p0014_main main.cpp)

# correctness: statement example + enumerated edges (n=1, path, star, before/after
# painting, repeated paints) + randomized differential against the BFS oracle.
add_executable(p0014_tests tests.cpp)
add_test(NAME p0014 COMMAND p0014_tests)

# worst case at n = m = 1e5: path (deepest centroid tree) and random tree.
add_executable(p0014_bench bench.cpp)
