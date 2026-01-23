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
include_directories(${OPRFSO_INSTALL_PATH}/include/secureJoin)

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
     ${CMAKE_SOURCE_DIR}/fuzzyPSI
)


set(COMMON_SOURCES ${SOURCES})

# handle build options for executables
option(BUILD_ALL "Build all executables" OFF)
option(BUILD_MIMT "Build mtest executable" OFF)
option(BUILD_OPRF  "Build oprf executable"  OFF)
option(BUILD_FMAP  "Build fmap executable"  ON)
option(BUILD_FPSI  "Build fpsi executable"  ON)

if(BUILD_ALL)
  set(BUILD_MIMT  ON)
  set(BUILD_OPRF   ON)
  set(BUILD_FMAP   ON)
  set(BUILD_FPSI   ON)
endif()

message(STATUS "BUILD_MIMT=${BUILD_MIMT}")
message(STATUS "BUILD_OPRF=${BUILD_OPRF}")
message(STATUS "BUILD_FMAP=${BUILD_FMAP}")
message(STATUS "BUILD_FPSI=${BUILD_FPSI}")

# create executables for tests and main program

if(BUILD_MIMT)
    add_executable(mtest
    ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/test.cpp
    ${COMMON_SOURCES}
    )
    target_include_directories(mtest PRIVATE ${HEADERS})
    target_link_libraries(mtest PRIVATE visa::volePSI visa::secureJoin)
endif()

if(BUILD_OPRF)
    add_executable(oprf
        ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/test_oprf.cpp
        ${COMMON_SOURCES}
    )
    target_include_directories(oprf PRIVATE ${HEADERS})
    target_link_libraries(oprf PRIVATE visa::volePSI visa::secureJoin)
endif()

if(BUILD_FMAP)
    add_executable(fmap
        ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/test_fmap.cpp
        ${COMMON_SOURCES}
    )
    target_include_directories(fmap PRIVATE ${HEADERS})
    target_link_libraries(fmap PRIVATE visa::volePSI visa::secureJoin)
endif()

if(BUILD_FPSI)
    add_executable(fpsi
        ${CMAKE_SOURCE_DIR}/fuzzyPSI/main/main.cpp
        ${COMMON_SOURCES}
    )
    target_include_directories(fpsi PRIVATE ${HEADERS})
    target_link_libraries(fpsi PRIVATE visa::volePSI visa::secureJoin)
endif()
