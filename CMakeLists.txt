cmake_minimum_required(VERSION 3.16)

project(libslater)

include(CTest)


set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

option(COVERAGE "Coverage" ON)
set(COMPILE_FLAGS -Wall -Werror -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -pthread)
set(LINK_FLAGS -fsanitize=address)


find_package(Boost 1.71 REQUIRED )


#SET(GCC_WERR_COMPILE_FLAGS "-Wall -Werror")
#SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${GCC_WERR_COMPILE_FLAGS} -fsanitize=address")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

find_package(PkgConfig)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${libslater_SOURCE_DIR} ${libslater_SOURCE_DIR}/include)

add_library(libslater SHARED src/libslater.cpp include/logger.h src/logger.cpp src/basis_functions.cpp src/engine.cpp include/homeier.h src/homeier.cpp src/parameters.cpp include/bfunctions.h src/bfunctions.cpp include/gaunt.h src/gaunt.cpp include/nested_summation.h src/nested_summation.cpp src/analytical-3c.cpp include/coordinates.h include/integrators.h include/analytical-3c.h)

target_compile_options(libslater PUBLIC ${COMPILE_FLAGS})
target_link_options(libslater PUBLIC ${LINK_FLAGS})

find_package(spdlog REQUIRED)
target_link_libraries(libslater PUBLIC spdlog::spdlog)
target_link_libraries (libslater PUBLIC Eigen3::Eigen)
target_link_libraries (libslater PUBLIC ${Boost_LIBRARIES})


add_executable(simple examples/simple.cpp)
add_dependencies(simple libslater)
target_include_directories(simple PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(simple PUBLIC ${COMPILE_FLAGS})
target_link_options(simple PUBLIC ${LINK_FLAGS})
target_link_libraries(simple libslater)

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
