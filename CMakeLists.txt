cmake_minimum_required (VERSION 3.25)

project (
    floxer
    LANGUAGES CXX
    VERSION 0.2.0
    DESCRIPTION "FM-Index longread PEX-seeded aligner"
)

# Force language standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ----- dependencies -----

# for the CPM package manager
# it is recommended to set the CPM_SOURCE_CACHE variable to avoid repetitive downloads
include(cmake/CPM.cmake)
CPMUsePackageLock (cmake/package-lock.cmake)

CPMGetPackage(sharg) # for CLI parsing
CPMGetPackage(spdlog) # for logging and output (includes libfmt)
CPMGetPackage(IVio) # for sequence input (that works well with fmindex-collection)
CPMGetPackage(IVSigma) # for sequence alphabet and transformation utilities (that works well with fmindex-collection)
CPMGetPackage(fmindex-collection) # for the FM-index
CPMGetPackage(cereal) # for FM-index serialization
CPMGetPackage(seqan3) # for alignment and sam/bam output
CPMGetPackage(interval-tree) # for interval tree data structure in the interval optimization
CPMGetPackage(use_ccache) # for faster compile times

# ----- floxer source code and build configuration -----

option (${PROJECT_NAME}_NO_OPENMP "Disable parallelism using OpenMP" OFF)
add_subdirectory(src/lib)
add_subdirectory(src/main)

# ----- testing -----

# for when the tests should not be built
option (${PROJECT_NAME}_TEST "Enable testing for ${PROJECT_NAME}." ON)

if (${PROJECT_NAME}_TEST)
    add_subdirectory (test EXCLUDE_FROM_ALL)
endif ()

# copy test data from source tree to build dir
set(test_data_dir "${CMAKE_BINARY_DIR}/test_data")
file (MAKE_DIRECTORY ${test_data_dir})

file (GLOB TEST_DATA_FILES ${PROJECT_SOURCE_DIR}/test/data/*)
foreach (test_data_file ${TEST_DATA_FILES})
    file(COPY ${test_data_file} DESTINATION ${test_data_dir})

    unset (test_data_file)
endforeach ()

unset(test_data_dir)
