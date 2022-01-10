cmake_minimum_required(VERSION 3.10)

project(wright VERSION 1.0
	DESCRIPTION "Computation of the Wright Function on the Real Line"
	HOMEPAGE_URL https://cirdans-home.github.io/mwright
LANGUAGES C Fortran)

configure_file(include/wright_config.h.in include/wright_config.h)

#------------------------------------------------------------------------------#
# Debug Options                                                                #
#------------------------------------------------------------------------------#
option(DEBUG "Set ON to enable debug prints" OFF)

# ---------------------------------------------------------------------------- #
# Compile the Library                                                          #
# ---------------------------------------------------------------------------- #
add_library(wright STATIC src/wright.c src/brentmod.f90 src/wrightmod.f90)
target_include_directories(wright PUBLIC
	"${PROJECT_BINARY_DIR}")
target_link_libraries(wright PUBLIC "-larb -lflint -lm")
if(DEBUG)
	add_definitions(-DDEBUG)
endif()

# ---------------------------------------------------------------------------- #
# Compile the Tests                                                            #
# ---------------------------------------------------------------------------- #
add_executable(optimalxitest test/optimalxitest.f90)
target_link_libraries(optimalxitest PUBLIC wright)

add_executable(wrighttest test/wrighttest.c)
target_link_libraries(wrighttest PUBLIC wright)
add_custom_command(
       TARGET wrighttest POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy
               ${CMAKE_SOURCE_DIR}/test/input.inp
						   ${CMAKE_CURRENT_BINARY_DIR}/input.inp)

add_executable(zwrighttest test/zwrighttest.c)
target_link_libraries(zwrighttest PUBLIC wright)

add_executable(mainarditest test/mainarditest.f90 src/brentmod.f90 src/wrightmod.f90)
target_link_libraries(wrighttest PUBLIC wright)

add_executable(quadtest test/quadtest.f90 src/brentmod.f90 src/wrightmod.f90)
target_link_libraries(wrighttest PUBLIC wright)

# ---------------------------------------------------------------------------- #
# Compile the docs                                                             #
# ---------------------------------------------------------------------------- #

option(doc_BUILD "Set to ON to build the documentation (needs Sphix)" ON)

if(doc_BUILD)
configure_file(docs/source/conf.py.in conf.py)
add_custom_command(
       TARGET wright PRE_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy
               ${CMAKE_BUILD_DIR}conf.py
               ${CMAKE_SOURCE_DIR}/docs/source/conf.py)

include(FindPackageHandleStandardArgs)

# We are likely to find Sphinx near the Python interpreter
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
    get_filename_component(_PYTHON_DIR "${PYTHON_EXECUTABLE}" DIRECTORY)
    set(
        _PYTHON_PATHS
        "${_PYTHON_DIR}"
        "${_PYTHON_DIR}/bin"
        "${_PYTHON_DIR}/Scripts")
endif()

find_program(
    SPHINX_EXECUTABLE
    NAMES sphinx-build sphinx-build.exe
    HINTS ${_PYTHON_PATHS})
mark_as_advanced(SPHINX_EXECUTABLE)

find_package_handle_standard_args(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)

if(Sphinx_FOUND)
    add_custom_command(
			TARGET wright POST_BUILD
			COMMAND make html
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs
			COMMENT "Building Sphinx documentation")
endif()
endif()
