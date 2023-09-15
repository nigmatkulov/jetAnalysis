cmake_minimum_required(VERSION 3.10)
project(simpleReadForest 
        VERSION 1.0 
        LANGUAGES C CXX)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED True)


message(STATUS "Looking for CERN ROOT")
#Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT CONFIG REQUIRED)
if (ROOT_FOUND)
        message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}") 
        message(STATUS "ROOT include dir found at ${ROOT_INCLUDE_DIRS}")
        message(STATUS "ROOT libraries found at ${ROOT_LIBRARIES}")
else()
        message(FATAL_ERROR "ROOT not found. CMake will exit!")

endif()

#Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# Directories where to search for include files
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
# Add same CXX flags that were used to combile and build ROOT (including CXX standard)
add_definitions(${ROOT_CXX_FLAGS})

add_executable(simpleReadForest simpleReadForest.C)
target_link_libraries(simpleReadForest ${ROOT_LIBRARIES})
target_include_directories(simpleReadForest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}) #${ROOT_INCLUDE_DIRS}

