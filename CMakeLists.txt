#-- require at least cmake version 3.8
cmake_minimum_required (VERSION 3.8)

#-- specify configuration options
option (USE_MPI "USE_MPI" OFF)
option (USE_OPENMP "USE_OPENMP" OFF)

#-- name the project and set the language
project (SuperLite)
set (SuperLite_VERSION_MAJOR 1)
set (SuperLite_VERSION_MINOR "x")
enable_language (Fortran)

#-- add the build directory to the include-file search path
include_directories(${PROJECT_BINARY_DIR})

#-- find MPI (cmake -DUSE_MPI=ON <src dir>)
find_package (MPI)
if (MPI_Fortran_FOUND AND USE_MPI)
  configure_file (mpimod_mpi.f mpimod.f COPYONLY)
else ()
  configure_file (mpimod_ser.f mpimod.f COPYONLY)
endif ()

#-- find OpenMP (cmake -DUSE_OPENMP=ON <src dir>)
find_package (OpenMP)
if (OPENMP_FOUND AND USE_OPENMP)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

#-- set the Fortran module directory
set (CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

#-- create lists of top-level source files
file (GLOB_RECURSE SuperLite_FILES "*.f" "*.f90")

#-- remove versions of mpimod from .f list
file (GLOB MPIMOD "mpimod_mpi.f" "mpimod_ser.f")
list (REMOVE_ITEM SuperLite_FILES ${MPIMOD})
list (APPEND SuperLite_FILES "${PROJECT_BINARY_DIR}/mpimod.f")

#-- add the executable and link MPI libraries
add_executable (superlite ${SuperLite_FILES})

#-- find LAPACK and BLAS ! added - gaw
find_package (BLAS)
find_package (LAPACK )
if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapckblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl")
  target_link_libraries(superlite PRIVATE ${lapackblas_libraries})
endif ()

#find_library (BLA_LIB libblas.so PATH /usr/lib/x86_64-linux-gnu/blas/)
#find_library (LAPACK_LIB liblapack.so PATH /usr/lib/x86_64-linux-gnu/lapack/)
#target_link_libraries(superlite PRIVATE ${BLA_LIB})
#target_link_libraries(superlite PRIVATE ${LAPACK_LIB})
#target_include_directories (superlite PUBLIC "/usr/lib/x86_64-linux-gnu/blas/")
#target_include_directories (superlite PUBLIC "/usr/lib/x86_64-linux-gnu/lapack/")

#target_link_libraries(superlite PRIVATE "/usr/lib/x86_64-linux-gnu/blas/libblas.so")
#target_link_libraries(superlite PRIVATE "/usr/lib/x86_64-linux-gnu/lapack/liblapack.so")
#target_include_directories (superlite PUBLIC "/usr/lib/x86_64-linux-gnu/blas/")
#target_include_directories (superlite PUBLIC "/usr/lib/x86_64-linux-gnu/lapack/")

#-- link library and include header for MPI
if (MPI_Fortran_FOUND AND USE_MPI)
  target_link_libraries (superlite PRIVATE ${MPI_C_LIBRARIES})
  target_include_directories (superlite PRIVATE ${MPI_C_INCLUDE_PATH})
  target_link_libraries (superlite PRIVATE ${MPI_Fortran_LIBRARIES})
  target_include_directories (superlite PRIVATE ${MPI_Fortran_INCLUDE_PATH})
endif ()

#-- temporary fix for MPI mismatch on latest MacOS
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")

#-- indicate successful configuration
message ("")
message ("SuperLite Version " ${SuperLite_VERSION_MAJOR} "."
	${SuperLite_VERSION_MINOR} " configured.")
message ("")
