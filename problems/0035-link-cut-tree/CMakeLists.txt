# the submittable driver.
add_executable(p0035_main main.cpp)

# correctness: hand shapes + edges + differential vs the brute path-walk oracle.
add_executable(p0035_tests tests.cpp)
add_test(NAME p0035 COMMAND p0035_tests)

# n=q=1e5 adversarial mix, measured.
add_executable(p0035_bench bench.cpp)
