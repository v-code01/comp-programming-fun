# the submittable driver: stdin -> A(x_1) … A(x_m).
add_executable(p0052_main main.cpp)

# correctness: hand examples + edges + the differential vs Horner on thousands
# of small instances (m<n, m>n, repeats, x=0).
add_executable(p0052_tests tests.cpp)
add_test(NAME p0052 COMMAND p0052_tests)

# product tree + remainder tree wall time at n=m=1e5.
add_executable(p0052_bench bench.cpp)
