cmake_minimum_required(VERSION 3.15)
project(MyProj)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)

# set install paths for volePSI and secure-join
set(VOLEPSI_INSTALL_PATH "${CMAKE_SOURCE_DIR}/install/volepsi")
set(OPRFSO_INSTALL_PATH "${CMAKE_SOURCE_DIR}/install/secure-join")

list(APPEND CMAKE_PREFIX_PATH
    "${VOLEPSI_INSTALL_PATH}"
    "${OPRFSO_INSTALL_PATH}"
)

# find volePSI and secure-join packages, and add include directories
find_package(volePSI REQUIRED) 
find_package(secureJoin REQUIRED)
include("${OPRFSO_INSTALL_PATH}/lib/cmake/secureJoin/secureJoinTargets.cmake")

include_directories(${VOLEPSI_INSTALL_PATH}/include)
include_directories(${OPRFSO_INSTALL_PATH}/include)

# Enable debug logging option
option(ENABLE_DEBUG_LOG "Enable debug logging" OFF)
if(ENABLE_DEBUG_LOG)
    add_compile_definitions(ENABLE_DEBUG_LOG=1)
    message(STATUS "Debug logging ENABLED")
else()
    add_compile_definitions(ENABLE_DEBUG_LOG=0)
    message(STATUS "Debug logging DISABLED")
endif()

# cmake build type and compiler options
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-O3)
endif()


# set source and header files
set(FUZZYPSI_DIR      ${CMAKE_SOURCE_DIR}/fuzzyPSI/psi)
set(MPEQT_DIR      ${CMAKE_SOURCE_DIR}/fuzzyPSI/mPeqt)
set(OPRFSO_DIR      ${CMAKE_SOURCE_DIR}/fuzzyPSI/oprf_so)

set(SOURCES
    ${FUZZYPSI_DIR}/psiLinfty.cpp
    ${FUZZYPSI_DIR}/psiL1.cpp
    ${FUZZYPSI_DIR}/psiL2.cpp
    ${FUZZYPSI_DIR}/mImt.cpp
    ${FUZZYPSI_DIR}/fmap.cpp
    ${MPEQT_DIR}/peqt.cpp
    ${OPRFSO_DIR}/oprfso.cpp
)

set(HEADERS
    ${FUZZYPSI_DIR}/psi.h
    ${FUZZYPSI_DIR}/mImt.h
    ${FUZZYPSI_DIR}/Defines.h
    ${FUZZYPSI_DIR}/fmap.h
    ${MPEQT_DIR}/peqt.h
    ${OPRFSO_DIR}/oprfso.h
    ${MPEQT_DIR}/n1not.h
)

set(DEBUG_HEADER ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/debug.h)

set(COMMON_SOURCES
    ${DEBUG_HEADER}
    ${SOURCES}
    ${HEADERS}
)


# create executables for tests and main program
add_executable(mtest
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/test.cpp
    ${COMMON_SOURCES}
)
target_include_directories(mtest PRIVATE
    ${CMAKE_SOURCE_DIR}/fuzzyPSI
    ${OPRFSO_INSTALL_PATH}/include/secureJoin
)
target_link_libraries(mtest PRIVATE visa::volePSI visa::secureJoin)

add_executable(m_oprf_test
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/test_oprf.cpp
    ${COMMON_SOURCES}
)
target_include_directories(m_oprf_test PRIVATE
    ${CMAKE_SOURCE_DIR}/fuzzyPSI
    ${OPRFSO_INSTALL_PATH}/include/secureJoin
)
target_link_libraries(m_oprf_test PRIVATE visa::volePSI visa::secureJoin)

add_executable(m_fmap_test
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/test_fmap.cpp
    ${COMMON_SOURCES}
)
target_include_directories(m_fmap_test PRIVATE
    ${CMAKE_SOURCE_DIR}/fuzzyPSI
    ${OPRFSO_INSTALL_PATH}/include/secureJoin
)
target_link_libraries(m_fmap_test PRIVATE visa::volePSI visa::secureJoin)

add_executable(main
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/main.cpp
    ${COMMON_SOURCES}
)
target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}/fuzzyPSI
    ${OPRFSO_INSTALL_PATH}/include/secureJoin  
)
target_link_libraries(main PRIVATE visa::volePSI visa::secureJoin)

set(ALL_TESTS mtest m_fmap_test)

foreach(target IN LISTS ALL_TESTS)
    set_target_properties(${target} PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
    )
    target_compile_options(${target} PRIVATE ${DEFAULT_COMPILER_OPTIONS_AND_WARNINGS})

    target_link_libraries(${target} PRIVATE
        visa::volePSI
        visa::secureJoin
    )

    add_test(NAME ${target} COMMAND ${target})
endforeach()
