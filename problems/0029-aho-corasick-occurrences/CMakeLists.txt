# the submittable driver -- aho-corasick build + single text pass over stdin.
add_executable(p0029_main main.cpp)

# correctness: statement example + edges + randomized differential against the
# naive per-pattern scan oracle.
add_executable(p0029_tests tests.cpp)
add_test(NAME p0029 COMMAND p0029_tests)

# build cost at sum|P|=1e5 and feed cost at |T|=1e5.
add_executable(p0029_bench bench.cpp)
