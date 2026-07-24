# the submittable driver.
add_executable(p0057_main main.cpp)

# correctness: hand shapes + edges + differential vs the naive best-first oracle.
add_executable(p0057_tests tests.cpp)
add_test(NAME p0057 COMMAND p0057_tests)

# n=1e5 m=5e5 k=1e5 worst case, measured.
add_executable(p0057_bench bench.cpp)
