# one place for the flags every problem inherits. change it here, changes
# everywhere -- no per-problem drift.

# the standard. no fallback, no extensions. if the compiler can't do c++23,
# it can't build this repo.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# warnings are errors. a warning is a bug that hasn't bitten yet.
add_compile_options(-Wall -Wextra -Werror)

# native tuning, but only the spelling this compiler actually accepts. apple
# clang wants -mcpu=native on arm; gcc/clang elsewhere want -march=native.
# guess wrong and the whole build dies on flag one -- so ask first.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" CPF_HAS_MARCH_NATIVE)
if(CPF_HAS_MARCH_NATIVE)
    set(CPF_TUNE -march=native)
else()
    check_cxx_compiler_flag("-mcpu=native" CPF_HAS_MCPU_NATIVE)
    if(CPF_HAS_MCPU_NATIVE)
        set(CPF_TUNE -mcpu=native)
    else()
        set(CPF_TUNE)
    endif()
endif()

# release is the real thing -- full opt, tuned for this machine. this is what
# the benchmarks measure.
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O2 ${CPF_TUNE})
endif()

# a build type for the things -Werror can't see until runtime: use-after-free,
# signed overflow, ub. tests live here. -O1 keeps it honest without hiding bugs
# behind aggressive inlining.
if(CMAKE_BUILD_TYPE STREQUAL "Sanitize")
    add_compile_options(-O1 -g -fno-omit-frame-pointer
                        -fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()
