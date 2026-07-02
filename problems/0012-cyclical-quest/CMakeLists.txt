# the submittable driver -- SAM build + doubled-walk queries over stdin.
add_executable(p0012_main main.cpp)

# correctness: statement example + edges + exhaustive small grid + randomized
# differential against the rotation-set oracle.
add_executable(p0012_tests tests.cpp)
add_test(NAME p0012 COMMAND p0012_tests)

# build cost at |s|=1e6 and query cost at sum|x|~1e6.
add_executable(p0012_bench bench.cpp)
