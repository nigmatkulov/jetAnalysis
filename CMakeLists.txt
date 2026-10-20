cmake_minimum_required(VERSION 3.10)
project(processForestSimple 
        VERSION 1.0 
        LANGUAGES C CXX)

option(DO_MIXED_PRIOR "Produce the optional mixed-prior unfolding histograms" OFF)

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
find_package(ROOT REQUIRED COMPONENTS RIO Net Physics Hist MathCore ROOTVecOps Rint Tree TreePlayer Core)
if (ROOT_FOUND)
        message(STATUS "ROOT ${ROOT_VERSION} found at ${ROOT_BINDIR}") 
        message(STATUS "ROOT include dir found at ${ROOT_INCLUDE_DIRS}")
        message(STATUS "ROOT libraries found at ${ROOT_LIBRARIES}")
else()
        message(FATAL_ERROR "ROOT not found. CMake will exit!")

endif()

include(${ROOT_USE_FILE})

# Directories where to search for include files
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
# Add same CXX flags that were used to combile and build ROOT (including CXX standard)
add_definitions(${ROOT_CXX_FLAGS})

# Generate ROOT dictionary for helper classes that use ClassDef/ClassImp
root_generate_dictionary(G__processForestSimple
                         JetCorrector.h
                         JetUncertainty.h
                         LINKDEF LinkDef.h)

add_library(jetanalysis_corrections STATIC
    JetCorrector.cc
    JetUncertainty.cc
    ${CMAKE_CURRENT_BINARY_DIR}/G__processForestSimple.cxx
)
set_target_properties(jetanalysis_corrections PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(jetanalysis_corrections PUBLIC ${ROOT_LIBRARIES})
target_include_directories(jetanalysis_corrections PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_executable(processForestSimple processForestSimple.C)
target_link_libraries(processForestSimple jetanalysis_corrections)
target_include_directories(processForestSimple PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS}) #${ROOT_INCLUDE_DIRS}
if(DO_MIXED_PRIOR)
    target_compile_definitions(processForestSimple PRIVATE DO_MIXED_PRIOR)
endif()

include(CTest)
add_library(jetanalysis_framework SHARED
    framework/src/CMSForestReader.cpp
    framework/src/AnalysisManager.cpp)
add_library(JetAnalysis::Framework ALIAS jetanalysis_framework)
set_target_properties(jetanalysis_framework PROPERTIES EXPORT_NAME Framework)
target_compile_features(jetanalysis_framework PUBLIC cxx_std_17)
target_include_directories(jetanalysis_framework PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/framework/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(jetanalysis_framework PUBLIC ROOT::TreePlayer ROOT::Tree ROOT::RIO ROOT::Hist)

add_library(jetanalysis_dijet SHARED
    framework/src/DiJetAnalysis.cpp
    framework/src/Histograms.cpp
    framework/src/Output.cpp)
add_library(JetAnalysis::DiJet ALIAS jetanalysis_dijet)
set_target_properties(jetanalysis_dijet PROPERTIES EXPORT_NAME DiJet)
target_link_libraries(jetanalysis_dijet PUBLIC jetanalysis_framework PRIVATE jetanalysis_corrections)
if(DO_MIXED_PRIOR)
    target_compile_definitions(jetanalysis_dijet PRIVATE DO_MIXED_PRIOR)
endif()
add_executable(jetAnalysis framework/src/main.cpp)
target_link_libraries(jetAnalysis PRIVATE JetAnalysis::DiJet)
add_executable(compare_root_outputs framework/tools/compare_root_outputs.cpp)
target_link_libraries(compare_root_outputs ${ROOT_LIBRARIES})
add_executable(validate_forest_reader framework/tools/validate_reader.cpp)
target_link_libraries(validate_forest_reader PRIVATE JetAnalysis::Framework jetanalysis_corrections)
if(DO_MIXED_PRIOR)
    target_compile_definitions(validate_forest_reader PRIVATE DO_MIXED_PRIOR)
endif()
add_executable(check_forest_schema framework/tools/check_forest_schema.cpp)
target_include_directories(check_forest_schema PRIVATE framework/include)
target_link_libraries(check_forest_schema ${ROOT_LIBRARIES})
if(BUILD_TESTING)
    add_executable(test_dijet_composition framework/tests/test_dijet_composition.cpp)
    target_link_libraries(test_dijet_composition PRIVATE JetAnalysis::DiJet)
    add_executable(test_analysis_manager framework/tests/test_analysis_manager.cpp)
    target_link_libraries(test_analysis_manager PRIVATE JetAnalysis::Framework)
    add_test(NAME analysis_manager COMMAND test_analysis_manager)
    find_program(JETANALYSIS_PYTHON NAMES python python3 HINTS "${CMAKE_CURRENT_SOURCE_DIR}/py-env/bin")
    if(JETANALYSIS_PYTHON)
        add_test(NAME library_consumer
            COMMAND ${JETANALYSIS_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/framework/tests/test_library_consumer.py
                    --build ${CMAKE_CURRENT_BINARY_DIR})
        add_test(NAME framework_workflow
            COMMAND ${JETANALYSIS_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/framework/tests/test_workflow.py
                    --build ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(framework_workflow PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    add_test(NAME framework_help COMMAND jetAnalysis --help)
    add_test(NAME framework_missing_arguments COMMAND jetAnalysis)
    set_tests_properties(framework_missing_arguments PROPERTIES WILL_FAIL TRUE)
    add_executable(test_root_branch_contract framework/tests/test_root_branch_contract.cpp)
    target_include_directories(test_root_branch_contract PRIVATE framework/include)
    target_link_libraries(test_root_branch_contract ${ROOT_LIBRARIES})
    add_test(NAME root_branch_contract COMMAND test_root_branch_contract)
    add_executable(test_analysis_contract tests/test_analysis_contract.cpp)
    target_include_directories(test_analysis_contract PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME analysis_contract COMMAND test_analysis_contract)
endif()

# A task application can consume installed libraries with find_package(), or use
# these same names from add_subdirectory() in a source checkout.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(framework/cmake/JetAnalysisFrameworkConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/JetAnalysisFrameworkConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/JetAnalysisFramework)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/JetAnalysisFrameworkConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(TARGETS jetanalysis_framework jetanalysis_dijet jetAnalysis
    EXPORT JetAnalysisFrameworkTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY framework/include/jetanalysis DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# The linked correction dictionary looks for this PCM beside the shared library.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libprocessForestSimple_rdict.pcm
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES JetCorrector.h JetUncertainty.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT JetAnalysisFrameworkTargets NAMESPACE JetAnalysis::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/JetAnalysisFramework)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/JetAnalysisFrameworkConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/JetAnalysisFrameworkConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/JetAnalysisFramework)
if(APPLE)
    set_target_properties(jetAnalysis PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    set_target_properties(jetanalysis_dijet jetanalysis_framework PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
    set_target_properties(jetAnalysis PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    set_target_properties(jetanalysis_dijet jetanalysis_framework PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()
