cmake_minimum_required(VERSION 3.16)

project(sttt_monitoring LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

function(add_sttt_experiment target source)
    add_executable("${target}" "${source}")
    target_include_directories(
        "${target}"
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
    if(MSVC)
        target_compile_options(
            "${target}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:/O2>"
        )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(
            "${target}" PRIVATE "$<$<NOT:$<CONFIG:Debug>>:-O3>"
        )
    endif()
endfunction()

set(STTT_ONLINE_FORMULAS
    ac
    ad
    since
    response
    response01
    response02
)

foreach(mode IN ITEMS incremental naive)
    foreach(formula IN LISTS STTT_ONLINE_FORMULAS)
        add_sttt_experiment(
            "sttt_online_${mode}_${formula}"
            "benchmarks/online/${mode}/main_eval_${formula}.cpp"
        )
    endforeach()
endforeach()

add_sttt_experiment(
    sttt_offline_approximate_random
    benchmarks/offline/approximate/random.cpp
)
add_sttt_experiment(
    sttt_offline_approximate_mutual_separation
    benchmarks/offline/approximate/mutual_separation.cpp
)
add_sttt_experiment(
    sttt_offline_approximate_water_tanks
    benchmarks/offline/approximate/water_tanks.cpp
)

include(CTest)
if(BUILD_TESTING)
    add_sttt_experiment(
        sttt_test_online_monitoring
        tests/test_online_monitoring.cpp
    )
    add_sttt_experiment(
        sttt_test_timed_until
        tests/test_timed_until.cpp
    )

    add_sttt_experiment(
        sttt_test_case_study_approximate
        tests/test_case_study_approximate.cpp
    )
    add_sttt_experiment(
        sttt_test_case_study_numeric
        tests/test_case_study_numeric.cpp
    )
    add_sttt_experiment(
        sttt_test_untimed_unary
        tests/test_untimed_unary.cpp
    )
    add_sttt_experiment(
        sttt_test_timed_future_optimized
        tests/test_timed_future_optimized.cpp
    )
    add_sttt_experiment(
        sttt_test_untimed_past_optimized
        tests/test_untimed_past_optimized.cpp
    )
    add_sttt_experiment(
        sttt_test_timed_past_optimized
        tests/test_timed_past_optimized.cpp
    )
    if(MSVC)
        target_compile_options(sttt_test_online_monitoring PRIVATE /UNDEBUG)
        target_compile_options(sttt_test_timed_until PRIVATE /UNDEBUG)
        target_compile_options(sttt_test_case_study_approximate PRIVATE /UNDEBUG)
        target_compile_options(sttt_test_case_study_numeric PRIVATE /UNDEBUG)
        target_compile_options(sttt_test_untimed_unary PRIVATE /UNDEBUG)
        target_compile_options(
            sttt_test_timed_future_optimized PRIVATE /UNDEBUG
        )
        target_compile_options(
            sttt_test_untimed_past_optimized PRIVATE /UNDEBUG
        )
        target_compile_options(
            sttt_test_timed_past_optimized PRIVATE /UNDEBUG
        )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(sttt_test_online_monitoring PRIVATE -UNDEBUG)
        target_compile_options(sttt_test_timed_until PRIVATE -UNDEBUG)
        target_compile_options(sttt_test_case_study_approximate PRIVATE -UNDEBUG)
        target_compile_options(sttt_test_case_study_numeric PRIVATE -UNDEBUG)
        target_compile_options(sttt_test_untimed_unary PRIVATE -UNDEBUG)
        target_compile_options(
            sttt_test_timed_future_optimized PRIVATE -UNDEBUG
        )
        target_compile_options(
            sttt_test_untimed_past_optimized PRIVATE -UNDEBUG
        )
        target_compile_options(
            sttt_test_timed_past_optimized PRIVATE -UNDEBUG
        )
    endif()

    add_test(NAME online_monitoring COMMAND sttt_test_online_monitoring)
    add_test(NAME timed_until COMMAND sttt_test_timed_until)
    add_test(
        NAME case_study_approximate
        COMMAND sttt_test_case_study_approximate
    )
    add_test(
        NAME case_study_numeric
        COMMAND sttt_test_case_study_numeric
    )
    add_test(NAME untimed_unary COMMAND sttt_test_untimed_unary)
    add_test(
        NAME timed_future_optimized
        COMMAND sttt_test_timed_future_optimized
    )
    add_test(
        NAME untimed_past_optimized
        COMMAND sttt_test_untimed_past_optimized
    )
    add_test(
        NAME timed_past_optimized
        COMMAND sttt_test_timed_past_optimized
    )
endif()
