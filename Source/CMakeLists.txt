
# Cmake setup with import of necessary libraries

# Standard cmake version requirement
cmake_minimum_required(VERSION 2.6)

# Inform user
MESSAGE("Configuring program.")

# Project (fortran must be enabled for OpenCMISS compatibility)
PROJECT( MultigridPoisson LANGUAGES C CXX Fortran )

# Add all the local source files from the folder "src" to the makefile
file(GLOB_RECURSE MULTIGRID_POISSON_SRC_LIBRARY
    "*.hpp"
    "*.cpp"
)

# Application source: could make this to select various files
# or just compile different .exe's which would be a little weird
file(GLOB_RECURSE MULTIGRID_POISSON_SRC_APPLICATION
    "*.hpp"
    "*.cpp"
)


#---- external packages ----#
# 1. MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


# Add the execubale using the source files in (src)
add_executable(MultigridPoisson ${MULTIGRID_POISSON_SRC_APPLICATION})

# Link the local and external libraries to the executable
target_link_libraries( 
	MultigridPoisson
	mpi
	mpi_cxx
	${MPI_LIBRARIES}
)


# Install the application binary
install(TARGETS MultigridPoisson
        RUNTIME DESTINATION bin
)

# Install the header files
INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  FILES_MATCHING PATTERN "*.hpp"
  )

# Install the header files
INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  FILES_MATCHING PATTERN "*.h"
  )


INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  FILES_MATCHING PATTERN "*.m"
  )

INSTALL ( DIRECTORY ${CMAKE_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  FILES_MATCHING PATTERN "*.txt"
  )













































