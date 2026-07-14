# the submittable driver.
add_executable(p0040_main main.cpp)

# correctness: hand examples + edges + exhaustive small alphabets + randomized
# differential + structured shapes, all against the O(n^2) substring-set oracle.
add_executable(p0040_tests tests.cpp)
add_test(NAME p0040 COMMAND p0040_tests)

# wall-clock of the ukkonen build at |s|=2e5, worst cases included.
add_executable(p0040_bench bench.cpp)
