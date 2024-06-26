cmake_minimum_required(VERSION 3.22)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Specifies the build type on single-configuration generators.")
endif()

if(DEFINED PROJECT_NAME)
  set(SUBPROJECT ON)
  set(INSTALL_SUBPROJECTS ON CACHE BOOL "Install subproject dependencies")
else()
  set(SUBPROJECT OFF)
  set(INSTALL_SUBPROJECTS ON CACHE BOOL "Install subproject dependencies" FORCE)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

file(STRINGS VERSION CGMF_VER)
set(cgmf_version ${CGMF_VER} CACHE STRING "CGMF Library Version")

project(cgmf VERSION ${CGMF_VER} LANGUAGES C CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)

include(CMakeDependentOption)
include(GNUInstallDirs)
include(FetchContent)

CMAKE_DEPENDENT_OPTION(cgmf.installation
                       "Install CGMF (default is ON)" ON
                       "NOT SUBPROJECT" ${INSTALL_SUBPROJECTS}
)

mark_as_advanced(cgmf.installation)

add_subdirectory(libcgmf)

option(cgmf.tests "Configure tests for CGMF" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(cgmf.tests)
  enable_testing()
endif()


if(cgmf.installation)
  install(DIRECTORY data/
          DESTINATION "${CMAKE_INSTALL_DATADIR}/data/cgmf-${CGMF_VER}"
          FILES_MATCHING PATTERN "*.dat*"
  )
endif()

FetchContent_Declare(
  osiris GIT_REPOSITORY git@github.com:beykyle/osiris.git
  GIT_TAG "origin/main"
)
message("Fetching beykyle/osiris")

FetchContent_Declare(
  json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)
message("Fetching nlohmann/json")

FetchContent_MakeAvailable(json osiris)


option(BUILD_PYCGMF "BUILD_PYCGMF" OFF)
if(BUILD_PYCGMF)
  find_package(Python REQUIRED  COMPONENTS Interpreter Development NumPy)
  find_package(pybind11 REQUIRED CONFIG)
  add_subdirectory(pyCGMF)
else()
  # build executable only in traditional standalone projects
  if(NOT SUBPROJECT)
    add_subdirectory(utils)
  endif()
endif()

