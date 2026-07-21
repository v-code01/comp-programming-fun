# the submittable driver.
add_executable(p0047_main main.cpp)

# correctness: hand cases + edges + the differential against the brute oracle.
add_executable(p0047_tests tests.cpp)
add_test(NAME p0047 COMMAND p0047_tests)

# wall time of the O(nm) solve at n=1000, m=1e5.
add_executable(p0047_bench bench.cpp)
