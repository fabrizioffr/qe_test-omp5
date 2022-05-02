include(CheckFortranCompilerFlag)

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 2022.0.0)
	message(FATAL_ERROR "Intel OneAPI is supported from v2022.0.0")
endif()

# set optimization specific flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xCORE-AVX512 -qopt-zmm-usage=high")

if(QE_ENABLE_OPENMP_OFFLOAD)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fiopenmp -fopenmp-targets=spir64")
    # Not used now: this should be used together with qe_enable_omp_offload function
    set(QE_OMP_OFFLOAD_COMPILE_OPTIONS)
    list(APPEND QE_OMP_OFFLOAD_COMPILE_OPTIONS -fiopenmp -fopenmp-targets=spir64)
    #
    set(QE_OMP_OFFLOAD_LINK_OPTIONS)
    list(APPEND QE_OMP_OFFLOAD_LINK_OPTIONS -qmkl -fsycl -lsycl -lOpenCL -lmkl_sycl -liomp5)
    message("   Checking ifx OpenMP Offload related compile options: ${QE_OMP_OFFLOAD_COMPILE_OPTIONS}")
    message("   Checking ifx OpenMP Offload related link options: ${QE_OMP_OFFLOAD_LINK_OPTIONS}")
    # Checking of linking options fails
    #set(CMAKE_REQUIRED_LINK_OPTIONS ${QE_OMP_OFFLOAD_COMPILE_OPTIONS})
    #check_fortran_compiler_flag("${QE_OMP_OFFLOAD_COMPILE_OPTIONS}" IFX_VALID_OPTIONS)
    #unset(CMAKE_REQUIRED_LINK_OPTIONS)
    #if(NOT IFX_VALID_OPTIONS)
    #	unset(IFX_VALID_OPTIONS CACHE)
    #	message(FATAL_ERROR "ifx related option check failed! "
    #						"Please check CMakeError.log for the exact error.")
    #endif()
endif()
#list(APPEND IFX_OPTIONS ${IFX_COMPILE_OPTIONS})
#message("   ifx related compile options: ${IFX_OPTIONS}")
#set(CMAKE_REQUIRED_LINK_OPTIONS ${IFX_OPTIONS})
#check_fortran_compiler_flag("${IFX_OPTIONS}" IFX_VALID_OPTIONS)
#unset(CMAKE_REQUIRED_LINK_OPTIONS)
#if(NOT IFX_VALID_OPTIONS)
#	unset(IFX_VALID_OPTIONS CACHE)
#	message(FATAL_ERROR "ifx related option check failed! "
#						"Please check CMakeError.log for the exact error.")
#endif()
#
#add_library(compilerCustomConfig INTERFACE)
#target_compile_options(compilerCustomConfig
#	INTERFACE
#		${IFX_COMPILE_OPTIONS})
#target_link_options(compilerCustomConfig
#	INTERFACE
#		${IFX_COMPILE_OPTIONS}
#		${IFX_LINK_OPTIONS})
