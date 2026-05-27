cmake_minimum_required(VERSION 3.10)
project(processForestSimple 
        VERSION 1.0 
        LANGUAGES C CXX)

        # Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build: Debug, Release, RelWithDebInfo" FORCE)
endif()

# Print the build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Define common warning flags
set(WARNING_FLAGS
    -Wall
    -Wextra
    -Wpedantic
)

# Build type specific flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -DDEBUG")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
    elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -DNDEBUG")
    endif() 

    # Apply common warning flags
    add_compile_options(${WARNING_FLAGS})
endif()

message(STATUS "Looking for CERN ROOT")
#Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Physics Hist MathCore ROOTVecOps Rint Tree Core)
if (ROOT_FOUND)
        message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}") 
        message(STATUS "ROOT include dir found at ${ROOT_INCLUDE_DIRS}")
        message(STATUS "ROOT libraries found at ${ROOT_LIBRARIES}")
else()
        message(FATAL_ERROR "ROOT not found. CMake will exit!")

endif()

include(${ROOT_USE_FILE})

# Directories where to search for include files
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
# Add same CXX flags that were used to combile and build ROOT (including CXX standard)
add_definitions(${ROOT_CXX_FLAGS})

# Generate ROOT dictionary for helper classes that use ClassDef/ClassImp
root_generate_dictionary(G__processForestSimple
                         JetCorrector.h
                         JetUncertainty.h
                         LINKDEF LinkDef.h)

add_executable(processForestSimple
    processForestSimple.C
    JetCorrector.cc
    JetUncertainty.cc
    ${CMAKE_CURRENT_BINARY_DIR}/G__processForestSimple.cxx
)
target_link_libraries(processForestSimple ${ROOT_LIBRARIES})
target_include_directories(processForestSimple PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}) #${ROOT_INCLUDE_DIRS}
