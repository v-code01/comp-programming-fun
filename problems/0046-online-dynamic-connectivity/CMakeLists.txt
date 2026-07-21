# the submittable driver.
add_executable(p0046_main main.cpp)

# correctness: the ETT layer alone, then hand shapes + edges + differential vs the
# brute recompute oracle.
add_executable(p0046_tests tests.cpp)
add_test(NAME p0046 COMMAND p0046_tests)

# n=1e5, q=2e5 adversarial churn, measured.
add_executable(p0046_bench bench.cpp)
