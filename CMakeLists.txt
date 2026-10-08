cmake_minimum_required(VERSION 3.15)
project(fpsi_cmp)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)

# Set the volePSI install path and add its include directories.
set(VOLEPSI_INSTALL_PATH "${CMAKE_SOURCE_DIR}/install/volepsi")

list(APPEND CMAKE_PREFIX_PATH
    "${VOLEPSI_INSTALL_PATH}"
)

find_package(volePSI REQUIRED)

include_directories(${VOLEPSI_INSTALL_PATH}/include)

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

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
endif()


# set source and header files
set(FUZZYPSI_DIR ${CMAKE_SOURCE_DIR}/fuzzyPSI/psi)
set(MPEQT_DIR ${CMAKE_SOURCE_DIR}/fuzzyPSI/mPeqt)

set(SOURCES
    ${FUZZYPSI_DIR}/psiLinfty.cpp
    ${FUZZYPSI_DIR}/psiL1.cpp
    ${FUZZYPSI_DIR}/psiL2.cpp
    ${FUZZYPSI_DIR}/mImt.cpp
    ${FUZZYPSI_DIR}/fmap.cpp
    ${MPEQT_DIR}/peqt.cpp
    ${MPEQT_DIR}/sspeqt_support.cpp
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/ole/ot_ole.cpp
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/andpair/triple.cpp
)

set(HEADERS
    ${CMAKE_SOURCE_DIR}/fuzzyPSI
)


set(COMMON_SOURCES ${SOURCES})

# create executables for tests and main program
add_executable(fpsi
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/main.cpp
    ${COMMON_SOURCES}
)
target_include_directories(fpsi PRIVATE ${HEADERS})
target_link_libraries(fpsi PRIVATE visa::volePSI)