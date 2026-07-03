# the submittable driver -- HLD over stdin, fast IO.
add_executable(p0018_main main.cpp)

# correctness: statement example + enumerated edges (n=1, a==b, single edge,
# path, star, toggle discipline) + randomized differential against the path-walk
# oracle.
add_executable(p0018_tests tests.cpp)
add_test(NAME p0018 COMMAND p0018_tests)

# worst case at n = m = 1e5: path (deepest chains, stack-safety proof) and random
# tree.
add_executable(p0018_bench bench.cpp)
