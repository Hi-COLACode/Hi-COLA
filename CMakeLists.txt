cmake_minimum_required(VERSION 3.18)
project(HiCOLA LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
option(HICOLA_BUILD_NATIVE "Build the COLASolver executable" OFF)
option(HICOLA_USE_MPI "Enable MPI and MPI FFTW in COLASolver" OFF)

if(NOT HICOLA_BUILD_NATIVE)
    return()
endif()

find_package(GSL REQUIRED)
find_package(Lua REQUIRED)
find_path(FFTW_INCLUDE_DIR fftw3.h REQUIRED)
find_library(FFTW_LIBRARY NAMES fftw3 libfftw3-3 REQUIRED)

set(FML_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/HiCOLA/FML")
add_executable(hi-cola-nbody
    "${FML_ROOT}/COLASolver/src/Main.cpp"
    "${FML_ROOT}/Global/Global.cpp"
    "${FML_ROOT}/Spline/Spline.cpp"
    "${FML_ROOT}/ODESolver/ODESolver.cpp"
    "${FML_ROOT}/GadgetUtils/GadgetUtils.cpp"
    "${FML_ROOT}/ParameterMap/ParameterMap.cpp"
    "${FML_ROOT}/FileUtils/FileUtils.cpp"
    "${FML_ROOT}/Math/Math.cpp"
)
target_include_directories(hi-cola-nbody PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/HiCOLA"
    "${FFTW_INCLUDE_DIR}" "${LUA_INCLUDE_DIR}"
)
target_compile_definitions(hi-cola-nbody PRIVATE
    NDIM_SPACE=3 USE_FFTW USE_GSL USE_LUA
)
target_link_libraries(hi-cola-nbody PRIVATE
    "${FFTW_LIBRARY}" GSL::gsl GSL::gslcblas "${LUA_LIBRARIES}"
)
if(HICOLA_USE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    find_library(FFTW_MPI_LIBRARY NAMES fftw3_mpi REQUIRED)
    target_compile_definitions(hi-cola-nbody PRIVATE USE_MPI)
    target_link_libraries(hi-cola-nbody PRIVATE MPI::MPI_CXX "${FFTW_MPI_LIBRARY}")
endif()
if(WIN32)
    target_compile_definitions(hi-cola-nbody PRIVATE _USE_MATH_DEFINES NOMINMAX)
    target_link_libraries(hi-cola-nbody PRIVATE psapi)
endif()
install(TARGETS hi-cola-nbody RUNTIME DESTINATION HiCOLA/bin)
