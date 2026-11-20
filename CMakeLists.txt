cmake_minimum_required(VERSION 3.18)
project(HiCOLA LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
option(HICOLA_USE_MPI "Enable MPI and MPI FFTW in COLASolver" OFF)
option(HICOLA_USE_HEALPIX "Enable HEALPix lightcone density maps in COLASolver" OFF)

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
if(HICOLA_USE_HEALPIX)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(HEALPIX REQUIRED IMPORTED_TARGET healpix_cxx)
    # healpix_cxx headers contain OpenMP pragmas and its flags enable OpenMP.
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_compile_definitions(hi-cola-nbody PRIVATE USE_HEALPIX)
    target_link_libraries(hi-cola-nbody PRIVATE PkgConfig::HEALPIX OpenMP::OpenMP_CXX)
endif()
include(GNUInstallDirs)
install(TARGETS hi-cola-nbody RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
