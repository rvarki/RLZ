include(FetchContent)

cmake_minimum_required (VERSION 3.15)
project (rlz)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

###############################################################################
# Set install directory if not specified
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "Installation directory")
################################################################################
# Specify build type
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}")
################################################################################
# Locating Git
find_package(Git)
if(GIT_FOUND)
    message(STATUS "git found: ${GIT_EXECUTABLE}")
else()
    message(FATAL_ERROR "git not found. Cloning of submodules will not work.")
endif()
################################################################################
# Version number of project
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch of project: ${GIT_BRANCH}")
message(STATUS "Git commit hash of project: ${GIT_COMMIT_HASH}")
################################################################################
# Compiler and linker flags
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS}")
#################################################################################
# Check dependencies
message(STATUS "The source directory of the project is ${CMAKE_CURRENT_LIST_DIR}")

# include_directories has been legacy since 2012 
# https://stackoverflow.com/questions/13703647/how-to-properly-add-include-directories-with-cmake
set(PROJECT_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
set(THIRDPARTY_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)

message(STATUS "The include directories of the project is ${INCLUDE_DIRECTORIES}")

# Log directory for external content builds
make_directory(${CMAKE_CURRENT_BINARY_DIR}/logs)
set(INSTALL_FILE  ${CMAKE_CURRENT_BINARY_DIR}/install.cmake)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    link_libraries(${OpenMP_CXX_LIBRARIES})
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

# DIVSUFSORT
find_library(DIVSUFSORT_LIB divsufsort HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT_LIB)
    message(STATUS "divsufsort library not found. Building as an external content")
else()
    message(STATUS "divsufsort library found at ${DIVSUFSORT_LIB}.")
endif()

# DIVSUFSORT64
find_library(DIVSUFSORT64_LIB divsufsort64 HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT64_LIB)
    message(STATUS "divsufsort64 library not found. Building as an external content")
else()
    message(STATUS "divsufsort64 library found at ${DIVSUFSORT64_LIB}.")
endif()

# SDSL
find_library(SDSL_LIB sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
find_path(SDSL_SRC sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/include PATHS ${CONDA_PREFIX}/include)
if(NOT SDSL_LIB)
    message(STATUS "sdsl library not found. Building as an external content")
    FetchContent_Declare(sdsl
            GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
            # GIT_TAG "v2.1.1"
    )

    if(NOT sdsl_POPULATED)
        FetchContent_Populate(sdsl)
        # Install in the current binary dir
        execute_process(COMMAND ./install.sh ${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)
        file(APPEND ${INSTALL_FILE} "MESSAGE(STATUS \"Installing sdsl.\")\n")
        file(APPEND ${INSTALL_FILE} "execute_process(COMMAND ./install.sh \${CMAKE_INSTALL_PREFIX}
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)\n")
        if(result)
            message(FATAL_ERROR "CMake step for sdsl-lite failed: ${result}")
        endif()
    endif()

    set(SDSL_SRC ${CMAKE_CURRENT_BINARY_DIR}/include CACHE PATH "Path to include dir" FORCE)
    set(DIVSUFSORT_SRC ${CMAKE_CURRENT_BINARY_DIR}/include CACHE PATH "Path to include dir" FORCE)
    set(SDSL_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a CACHE PATH "Path to libsdsl.a" FORCE)
    set(DIVSUFSORT_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort.a CACHE PATH "Path to libdivsufsort.a" FORCE)
    set(DIVSUFSORT64_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort64.a CACHE PATH "Path to libdivsufsort64.a" FORCE)

else()
    message(STATUS "sdsl library found at ${SDSL_LIB}.")
    message(STATUS "sdsl sources found at ${SDSL_SRC}.")
endif()

add_library(fm_wrapper ${PROJECT_SRC_DIR}/fm_wrapper.cpp)
target_include_directories(fm_wrapper PUBLIC ${SDSL_SRC} ${DIVSUFSORT_SRC} ${DIVSUFSORT64_SRC} ${PROJECT_INCLUDE_DIR})
target_link_libraries(fm_wrapper PUBLIC ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB})

add_library(rlz_algo ${PROJECT_SRC_DIR}/rlz_algo.cpp)
target_include_directories(rlz_algo PUBLIC ${SDSL_SRC} ${DIVSUFSORT_SRC} ${DIVSUFSORT64_SRC} ${PROJECT_INCLUDE_DIR} ${THIRDPARTY_INCLUDE_DIR})
target_link_libraries(rlz_algo PUBLIC ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB} fm_wrapper)

add_library(rlz_algo_char ${PROJECT_SRC_DIR}/rlz_algo_char.cpp)
target_include_directories(rlz_algo_char PUBLIC ${SDSL_SRC} ${DIVSUFSORT_SRC} ${DIVSUFSORT64_SRC} ${PROJECT_INCLUDE_DIR} ${THIRDPARTY_INCLUDE_DIR})
target_link_libraries(rlz_algo_char PUBLIC ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB} fm_wrapper)

add_executable(rlz rlz.cpp)
target_include_directories(rlz PUBLIC ${SDSL_SRC} ${THIRDPARTY_INCLUDE_DIR} ${THIRDPARTY_INCLUDE_DIR}/CLI)
target_link_libraries(rlz PUBLIC ${SDSL_LIB} fm_wrapper rlz_algo rlz_algo_char)

install(TARGETS fm_wrapper rlz_algo rlz_algo_char TYPE LIBRARY)
install(TARGETS rlz TYPE RUNTIME)

