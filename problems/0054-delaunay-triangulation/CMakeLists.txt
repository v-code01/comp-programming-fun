# the submittable driver.
add_executable(p0054_main main.cpp)

# correctness: hand cases + the general-position gate + the differential against
# the O(n^4) empty-circle Delaunay-edge oracle.
add_executable(p0054_tests tests.cpp)
add_test(NAME p0054 COMMAND p0054_tests)

# throughput of Bowyer-Watson at n=2000.
add_executable(p0054_bench bench.cpp)
