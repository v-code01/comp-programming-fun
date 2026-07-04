# the submittable driver -- read s and t, print the LCS length.
add_executable(p0028_main main.cpp)

# correctness: hand examples + seam edges + exhaustive small grid + randomized
# differential against the O(nm) LCS-DP oracle.
add_executable(p0028_tests tests.cpp)
add_test(NAME p0028 COMMAND p0028_tests)

# suffix-array build cost at |w| = 2e5, worst-case and random shapes.
add_executable(p0028_bench bench.cpp)
