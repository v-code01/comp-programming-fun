# the submittable driver: stdin -> B = A^{-1} mod x^n.
add_executable(p0043_main main.cpp)

# correctness: hand examples + edges + the differential vs the direct O(n^2)
# inverse AND the A*B==1 verifier, on thousands of small polynomials.
add_executable(p0043_tests tests.cpp)
add_test(NAME p0043 COMMAND p0043_tests)

# newton + NTT wall time at n=1e5.
add_executable(p0043_bench bench.cpp)
