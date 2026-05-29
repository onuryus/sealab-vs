cmake_minimum_required(VERSION 3.16)
project(sealab-vs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# ---------------------------------------------------------------------------
# If a conda env is active, prefer its compiler/Boost/RDKit (ABI consistency).
# Put CONDA_PREFIX at the front of the search path.
# ---------------------------------------------------------------------------
if(DEFINED ENV{CONDA_PREFIX})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
  set(CONDA_LIB "$ENV{CONDA_PREFIX}/lib")
  set(CONDA_RDKIT_INC "$ENV{CONDA_PREFIX}/include/rdkit")
  message(STATUS "conda env detected: $ENV{CONDA_PREFIX}")
endif()

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)
find_package(Threads REQUIRED)

# Find RDKit headers explicitly: one predictable path that works for both
# conda and apt builds, independent of which exported target names the
# config package uses.
find_path(RDKIT_INCLUDE_DIR GraphMol/GraphMol.h
          HINTS ${CONDA_RDKIT_INC}
          PATHS /usr/include/rdkit /usr/local/include/rdkit)
if(NOT RDKIT_INCLUDE_DIR)
  message(FATAL_ERROR
    "RDKit headers not found.\n"
    "  conda:  conda install -c conda-forge librdkit-dev\n"
    "  apt:    sudo apt install librdkit-dev")
endif()
message(STATUS "RDKit headers: ${RDKIT_INCLUDE_DIR}")

set(RDKIT_COMPONENTS
    SmilesParse GraphMol Descriptors Fingerprints DataStructs
    RDGeneral RDGeometryLib SubstructMatch FilterCatalog Catalogs
    MolStandardize ChemTransforms FileParsers Subgraphs)

set(RDKIT_LIBS "")
foreach(comp ${RDKIT_COMPONENTS})
  find_library(RDKIT_${comp}_LIB NAMES RDKit${comp}
               HINTS ${CONDA_LIB}
               PATHS /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib)
  if(RDKIT_${comp}_LIB)
    list(APPEND RDKIT_LIBS ${RDKIT_${comp}_LIB})
  else()
    message(WARNING "RDKit component not found: RDKit${comp}")
  endif()
endforeach()

# ---------------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------------
add_executable(sealab-vs src/main.cpp src/chem.cpp)

target_include_directories(sealab-vs PRIVATE
    ${RDKIT_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
target_link_libraries(sealab-vs PRIVATE
    ${RDKIT_LIBS} OpenMP::OpenMP_CXX Threads::Threads)

# Embed the conda lib path as rpath so the binary finds its shared
# libraries even when the conda env is not activated (e.g. when launched
# from a GUI subprocess).
if(DEFINED CONDA_LIB)
  set_target_properties(sealab-vs PROPERTIES
      BUILD_RPATH "${CONDA_LIB}"
      INSTALL_RPATH "${CONDA_LIB}")
endif()
