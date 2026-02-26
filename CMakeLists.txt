cmake_minimum_required(VERSION 3.15)
project(ehijing LANGUAGES CXX)

# If you have FindPythia8.cmake in ./cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# ---- C++ standard (target-based, no global flags) ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Dependencies ----
find_package(Pythia8 REQUIRED)
find_package(HepMC3 REQUIRED)
find_package(GSL REQUIRED)

# ---- Executable ----
add_executable(ehijing
  src/main.cpp
)

# Includes (target-scoped)
target_include_directories(ehijing PRIVATE
  ${PYTHIA8_INCLUDE_DIR}
  ${HEPMC3_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
)

# Link libraries (target-scoped)
target_link_libraries(ehijing PRIVATE
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