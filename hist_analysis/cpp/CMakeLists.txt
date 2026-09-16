cmake_minimum_required(VERSION 3.10)

project(jetanalysis_regularization_scan
        VERSION 0.1.0
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the build type: Debug, Release, or RelWithDebInfo" FORCE)
endif()

find_package(Threads REQUIRED)

add_library(jetanalysis_regularization_scan SHARED
  regularization_scan.cpp
)

target_compile_features(jetanalysis_regularization_scan PUBLIC cxx_std_17)
target_include_directories(jetanalysis_regularization_scan
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(jetanalysis_regularization_scan PRIVATE Threads::Threads)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(jetanalysis_regularization_scan PRIVATE
    -Wall
    -Wextra
    -Wpedantic
  )
endif()

set_target_properties(jetanalysis_regularization_scan PROPERTIES
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME jetanalysis_regularization_scan
)

include(GNUInstallDirs)
install(TARGETS jetanalysis_regularization_scan
  EXPORT jetanalysis_regularization_scan_targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES regularization_scan.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
