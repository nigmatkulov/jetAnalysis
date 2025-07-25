cmake_minimum_required(VERSION 3.10)
project(JetAnalysis
        VERSION 1.0 
        LANGUAGES C CXX)

set(libname JetAnalysis)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

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


# find_package(ROOT CONFIG REQUIRED)
# if (ROOT_FOUND)
#         message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}") 
#         message(STATUS "ROOT include dir found at ${ROOT_INCLUDE_DIRS}")
#         message(STATUS "ROOT libraries found at ${ROOT_LIBRARIES}")
# else()
#         message(FATAL_ERROR "ROOT not found. CMake will exit!")
# endif()

#Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# Directories where to search for include files
include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR} 
        ${ROOT_INCLUDE_DIRS}
)

# Add same CXX flags that were used to combile and build ROOT (including CXX standard)
add_definitions(${ROOT_CXX_FLAGS})

# List all headers to be used to generate dictionary and added to the library
set(HEADERS
        JetESRAnalysis.h
        Event.h
        Track.h
        GenTrack.h
        RecoJet.h
        GenJet.h
        BaseAnalysis.h
        BaseReader.h
        BaseHistoManager.h
        BaseJet.h
        DiJet.h
        HistoManagerJetESR.h
        HistoManagerDiJet.h
        Collections.h
        ForestAODReader.h
        JetCorrector.h
        JetUncertainty.h
        EventCut.h
        JetCut.h
        DiJetCut.h
        Manager.h
        TriggerAndSkim.h
        DiJetAnalysis.h
)

# List source files
set(SRC 
        JetESRAnalysis.cc
        Event.cc
        Track.cc
        GenTrack.cc
        RecoJet.cc
        GenJet.cc
        BaseReader.cc
        BaseJet.cc
        DiJet.cc
        HistoManagerJetESR.cc
        HistoManagerDiJet.cc
        ForestAODReader.cc
        JetCorrector.cc
        JetUncertainty.cc
        EventCut.cc
        JetCut.cc
        DiJetCut.cc
        Manager.cc
        TriggerAndSkim.cc
        DiJetAnalysis.cc
)

# Generate ROOT dictionaries
root_generate_dictionary(G__${libname} 
                         ${HEADERS}
                         LINKDEF LinkDef.h
)

# Generate rootmap
root_generate_rootmap(${libname} LINKDEF LinkDef.h)

#ROOT_LINKER_LIBRARY(${libname} ${SRC} G__${libname}.cxx)

# Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SRC} G__${libname}.cxx)
# Link libraries need to run the code
target_link_libraries(${libname} ${ROOT_LIBRARIES})

# install(TARGETS $(libname)
#         LIBRARY DESTINATION "${CMAKE_BINARY_DIR}/lib" OPTIONAL)

# Create main executable (to call from the command line)
add_executable(jetESR jetESR.cxx)
# Link created libraries
target_link_libraries(jetESR ${libname})

# Create main executable (to call from the command line)
add_executable(dijetAna_pPb8160 dijetAna_pPb8160.cxx)
# Link created libraries
target_link_libraries(dijetAna_pPb8160 ${libname})

# Create main executable (to call from the command line)
add_executable(dijetAna_pPb5020 dijetAna_pPb5020.cxx)
# Link created libraries
target_link_libraries(dijetAna_pPb5020 ${libname})

# Create main executable (to call from the command line)
add_executable(dijetAna_pp5020 dijetAna_pp5020.cxx)
# Link created libraries
target_link_libraries(dijetAna_pp5020 ${libname})

# Include directories 
#target_include_directories(jetAna PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

# (Optional) Print final compiler flags
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
