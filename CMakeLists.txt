cmake_minimum_required(VERSION 3.15)
project(ehijing LANGUAGES CXX)

option(EHIJING_BUILD_DOCS_ONLY "Configure only the eHIJING Doxygen documentation target" OFF)

# ---- Documentation ----
find_package(Doxygen QUIET)

if(DOXYGEN_FOUND)
  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating eHIJING Doxygen documentation"
    VERBATIM
  )
else()
  add_custom_target(docs
    COMMAND ${CMAKE_COMMAND} -E echo "Doxygen was not found. Install doxygen, then rerun: cmake --build <build-dir> --target docs"
    COMMAND ${CMAKE_COMMAND} -E false
    COMMENT "Doxygen is required to generate eHIJING documentation"
    VERBATIM
  )
endif()

if(EHIJING_BUILD_DOCS_ONLY)
  return()
endif()

# If you have FindPythia8.cmake in ./cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# ---- Dependencies ----
find_package(Pythia8 REQUIRED)
find_package(HepMC3 REQUIRED)
find_package(GSL REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
if(TARGET yaml-cpp::yaml-cpp)
  set(EHIJING_YAML_TARGET yaml-cpp::yaml-cpp)
else()
  set(EHIJING_YAML_TARGET yaml-cpp)
endif()

include(GNUInstallDirs)
set(EHIJING_NUCLEUS_DATABASE "${CMAKE_CURRENT_SOURCE_DIR}/../../input/ehijing/woods_saxon.yaml"
    CACHE FILEPATH "Nuclear parameter YAML to use and install")
add_library(ehijing_config STATIC src/cli.cpp src/nucleus_database.cpp)
target_compile_features(ehijing_config PUBLIC cxx_std_17)
target_include_directories(ehijing_config PUBLIC include)
target_link_libraries(ehijing_config PUBLIC ${EHIJING_YAML_TARGET})
target_compile_definitions(ehijing_config PRIVATE
  EHIJING_SOURCE_NUCLEUS_DATABASE="${EHIJING_NUCLEUS_DATABASE}"
  EHIJING_INSTALLED_NUCLEUS_DATABASE="${CMAKE_INSTALL_FULL_DATADIR}/electra/ehijing/woods_saxon.yaml"
)

# ---- Executable ----
add_executable(ehijing
  src/file_digest.cpp
  src/dis_kinematics.cpp
  src/event_output.cpp
  src/filesystem_utils.cpp
  src/hadronizer.cpp
  src/main.cpp
  src/modified_ff.cpp
  src/spacetime_history.cpp
)

# ---- C++ standard (target-based) ----
target_compile_features(ehijing PRIVATE cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Includes (target-scoped)
target_include_directories(ehijing PRIVATE
  include
  ${PYTHIA8_INCLUDE_DIR}
  ${HEPMC3_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
)

# Link libraries (target-scoped)
target_link_libraries(ehijing PRIVATE
  ehijing_config
  ${PYTHIA8_LIBRARIES}
  ${HEPMC3_LIBRARIES}
  ${GSL_LIBRARIES}
)

# Threads: portable replacement for -pthread
find_package(Threads REQUIRED)
target_link_libraries(ehijing PRIVATE Threads::Threads)

# You usually don't need to link -ldl manually on macOS; keep it Linux-only
if(UNIX AND NOT APPLE)
  target_link_libraries(ehijing PRIVATE dl)
endif()

# Install
include(GNUInstallDirs)
install(TARGETS ehijing RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(EXISTS "${EHIJING_NUCLEUS_DATABASE}")
  install(FILES "${EHIJING_NUCLEUS_DATABASE}"
          DESTINATION ${CMAKE_INSTALL_DATADIR}/electra/ehijing RENAME woods_saxon.yaml)
endif()

# Synthetic diagnostics tests: no event generation or reference data required.
include(CTest)
if(BUILD_TESTING)
  add_executable(nucleus_database_test tests/nucleus_database_test.cpp)
  target_link_libraries(nucleus_database_test PRIVATE ehijing_config)
  add_test(NAME nucleus_database COMMAND nucleus_database_test --test-database "${EHIJING_NUCLEUS_DATABASE}")
  add_executable(seed_policy_test tests/seed_policy_test.cpp src/file_digest.cpp)
  target_compile_features(seed_policy_test PRIVATE cxx_std_17)
  target_include_directories(seed_policy_test PRIVATE include ${PYTHIA8_INCLUDE_DIR})
  target_link_libraries(seed_policy_test PRIVATE ehijing_config ${PYTHIA8_LIBRARIES} ${GSL_LIBRARIES} Threads::Threads)
  add_test(NAME seed_policy COMMAND seed_policy_test)
  add_executable(spacetime_history_test
    tests/spacetime_history_test.cpp
    src/spacetime_history.cpp
    src/filesystem_utils.cpp
    src/event_output.cpp
  )
  target_compile_features(spacetime_history_test PRIVATE cxx_std_17)
  target_include_directories(spacetime_history_test PRIVATE
    include ${PYTHIA8_INCLUDE_DIR} ${HEPMC3_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
  )
  target_link_libraries(spacetime_history_test PRIVATE
    ${PYTHIA8_LIBRARIES} ${HEPMC3_LIBRARIES} ${GSL_LIBRARIES} Threads::Threads
  )
  add_test(NAME spacetime_history COMMAND spacetime_history_test)
  add_executable(spectator_balance_test
    tests/spectator_balance_test.cpp
    src/spacetime_history.cpp
    src/event_output.cpp
  )
  target_compile_features(spectator_balance_test PRIVATE cxx_std_17)
  target_include_directories(spectator_balance_test PRIVATE
    include ${PYTHIA8_INCLUDE_DIR} ${HEPMC3_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
  )
  target_link_libraries(spectator_balance_test PRIVATE
    ${PYTHIA8_LIBRARIES} ${HEPMC3_LIBRARIES} ${GSL_LIBRARIES} Threads::Threads
  )
  add_test(NAME spectator_balance COMMAND spectator_balance_test)
  add_executable(nuclear_target_test tests/nuclear_target_test.cpp)
  target_compile_features(nuclear_target_test PRIVATE cxx_std_17)
  target_include_directories(nuclear_target_test PRIVATE include ${PYTHIA8_INCLUDE_DIR})
  target_link_libraries(nuclear_target_test PRIVATE ${PYTHIA8_LIBRARIES} ${GSL_LIBRARIES} Threads::Threads)
  add_test(NAME nuclear_target COMMAND nuclear_target_test)
  add_executable(shared_remnant_test tests/shared_remnant_test.cpp src/modified_ff.cpp src/spacetime_history.cpp)
  target_compile_features(shared_remnant_test PRIVATE cxx_std_17)
  target_include_directories(shared_remnant_test PRIVATE include ${PYTHIA8_INCLUDE_DIR} ${HEPMC3_INCLUDE_DIR} ${GSL_INCLUDE_DIR})
  target_link_libraries(shared_remnant_test PRIVATE ${PYTHIA8_LIBRARIES} ${HEPMC3_LIBRARIES} ${GSL_LIBRARIES} Threads::Threads)
  add_test(NAME shared_remnant COMMAND shared_remnant_test)
endif()
