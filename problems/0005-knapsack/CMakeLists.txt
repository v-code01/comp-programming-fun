# the submittable driver.
add_executable(p0005_main main.cpp)

# correctness: statement examples, edges, and differentials against two oracles.
add_executable(p0005_tests tests.cpp)
add_test(NAME p0005 COMMAND p0005_tests)

# cost of the take-all shortcut and the large-branch dp.
add_executable(p0005_bench bench.cpp)
