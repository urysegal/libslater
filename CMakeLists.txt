cmake_minimum_required(VERSION 3.16)

project(libslater)

include(CTest)


set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

option(COVERAGE "Coverage" ON)
#set(COMPILE_FLAGS -Wall -Werror -fno-omit-frame-pointer -pthread)

set(COMPILE_FLAGS -Wall -Werror -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -pthread)
set(LINK_FLAGS -fsanitize=address)


find_package(Boost 1.71 REQUIRED )


find_package(PkgConfig)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${libslater_SOURCE_DIR} ${libslater_SOURCE_DIR}/include)

add_library(libslater SHARED src/libslater.cpp include/logger.h src/logger.cpp src/basis_functions.cpp src/engine.cpp
        include/homeier.h src/homeier.cpp src/parameters.cpp include/bfunctions.h src/bfunctions.cpp include/gaunt.h
        src/gaunt.cpp include/nested_summation.h src/nested_summation.cpp src/analytical-3c.cpp include/coordinates.h
        include/integrators.h include/analytical-3c.h precalculate/gaunt-table.h precalculate/gaunt-table.cpp
        include/slater-utils.h src/utils.cpp src/semi-infinite-3c.cpp include/nested_summation_debug.h

        src/fortran-4c/bessel_reduite.c
        src/fortran-4c/binom.c
        src/fortran-4c/calf42.c
        src/fortran-4c/cartco.c
        src/fortran-4c/cnp.c
        src/fortran-4c/coeffbons.c
        src/fortran-4c/constlmp.c
        src/fortran-4c/Dbsj.c
        src/fortran-4c/dfactorial.c
        src/fortran-4c/DFj42.c
        src/fortran-4c/DG.c
        src/fortran-4c/dobin.c
        src/fortran-4c/Dtk.c
        src/fortran-4c/DX.c
        src/fortran-4c/fa42.c
        src/fortran-4c/factorial.c
        src/fortran-4c/Fj42.c
        src/fortran-4c/FSCj.c
        src/fortran-4c/gaunt1.c
        src/fortran-4c/gaunt.c
        src/fortran-4c/GaussLag.c
        src/fortran-4c/GaussLeg.c
        src/fortran-4c/g.c
        src/fortran-4c/GREP-WD42.c
        src/fortran-4c/harmonique.c
        src/fortran-4c/hatk.c
        src/fortran-4c/HFa42.c
        src/fortran-4c/HFj42.c
        src/fortran-4c/legendre.c
        src/fortran-4c/LU.c
        src/fortran-4c/main.c
        src/fortran-4c/miscfunc.c
        src/fortran-4c/msphj.c
        src/fortran-4c/power_x.c
        src/fortran-4c/racdi.c
        src/fortran-4c/sdbar42_lrzv.c
        src/fortran-4c/sinf42.c
        src/fortran-4c/vcoord.c
        src/fortran-4c/vharmonique.c
        src/fortran-4c/zero_bsj.c

        include/safouhi-4c.h src/safouhi-4c.cpp)

target_compile_options(libslater PUBLIC ${COMPILE_FLAGS})
target_link_options(libslater PUBLIC ${LINK_FLAGS})

find_package(spdlog REQUIRED)
target_link_libraries(libslater PUBLIC spdlog::spdlog)
target_link_libraries (libslater PUBLIC Eigen3::Eigen)
target_link_libraries (libslater PUBLIC ${Boost_LIBRARIES})
target_link_libraries (libslater PUBLIC f2c m)


add_executable(simple examples/simple.cpp)
add_dependencies(simple libslater)
target_include_directories(simple PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(simple PUBLIC ${COMPILE_FLAGS})
target_link_options(simple PUBLIC ${LINK_FLAGS})
target_link_libraries(simple libslater)

add_executable(test3c examples/test3c.cpp)
add_dependencies(test3c libslater)
target_include_directories(test3c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(test3c PUBLIC ${COMPILE_FLAGS})
target_link_options(test3c PUBLIC ${LINK_FLAGS})
target_link_libraries(test3c libslater)


if(COVERAGE)
    message(STATUS "Enabling coverage")
    target_compile_options(libslater PUBLIC --coverage)
    target_link_options(libslater PUBLIC --coverage)

    find_program(CODECOV_GCOV gcov)
    find_program(CODECOV_LCOV lcov)
    find_program(CODECOV_GENHTML genhtml)

    set(CODECOV_OUTPUTFILE cmake_coverage.output)
    set(CODECOV_HTMLOUTPUTDIR coverage_results)

    add_custom_target(coverage_init ALL ${CODECOV_LCOV}  --base-directory .  --directory ${CMAKE_BINARY_DIR}/CMakeFiles/libslater.dir/src  --exclude '/usr/include/*' --output-file ${CODECOV_OUTPUTFILE} --capture --initial )
    add_custom_target(coverage ${CODECOV_LCOV}  --base-directory .  --directory ${CMAKE_BINARY_DIR}/CMakeFiles/libslater.dir/src --exclude '/usr/include/*' --output-file ${CODECOV_OUTPUTFILE} --capture COMMAND genhtml -o ${CODECOV_HTMLOUTPUTDIR} ${CODECOV_OUTPUTFILE} )
endif()

add_subdirectory(test)

add_custom_target(cleancov
        COMMAND find ${CMAKE_BINARY_DIR} -name '*.gc??' -delete
        )

set_target_properties(libslater PROPERTIES PUBLIC_HEADER libslater.h)
include(GNUInstallDirs)
install(TARGETS libslater
        PUBLIC_HEADER
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libslater
        )



option(BUILD_DOC "Build documentation" )

if (BUILD_DOC)
# check if Doxygen is installed
find_package(Doxygen REQUIRED dot)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
endif()
