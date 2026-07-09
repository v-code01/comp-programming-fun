# the submittable driver.
add_executable(p0039_main main.cpp)

# correctness: hand shapes + edges + differential vs the brute recompute oracle.
add_executable(p0039_tests tests.cpp)
add_test(NAME p0039 COMMAND p0039_tests)

# n=q=2e5 adversarial add/remove churn, measured.
add_executable(p0039_bench bench.cpp)
