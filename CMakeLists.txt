cmake_minimum_required (VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (svgfill)

if(WASM_BUILD)
    set(CMAKE_FIND_ROOT_PATH_BACKUP "${CMAKE_FIND_ROOT_PATH}")
    set(CMAKE_FIND_ROOT_PATH "")
endif()

find_package(Boost REQUIRED CONFIG)
find_package(svgpp REQUIRED CONFIG)
find_package(libxml2 REQUIRED CONFIG)
find_package(CGAL REQUIRED CONFIG)
find_package(gmp REQUIRED CONFIG)
find_package(mpfr REQUIRED CONFIG)


if(WIN32 AND ("$ENV{CONDA_BUILD}" STREQUAL ""))
	set(Boost_USE_STATIC_LIBS ON)
	set(Boost_USE_MULTITHREADED ON)
    if (USE_STATIC_MSVC_RUNTIME)
        set(Boost_USE_STATIC_RUNTIME ON)
    endif()
else()
    # Disable Boost's autolinking as the libraries to be linked to are supplied
    # already by CMake, and it's going to conflict if there are multiple, as is
    # the case in conda-forge's libboost feedstock.
    add_definitions(-DBOOST_ALL_NO_LIB)
    if(WIN32)
        # Necessary for boost version >= 1.67
        set(BCRYPT_LIBRARIES "bcrypt.lib")
    endif()
endif()

if (MSVC)
	add_definitions(-bigobj)
endif()

if(WASM_BUILD)
    set(CMAKE_FIND_ROOT_PATH "${CMAKE_FIND_ROOT_PATH_BACKUP}")
endif()


file(GLOB LIB_H_FILES src/*.h)
message(STATUS "Found header files: ${LIB_H_FILES}")
file(GLOB LIB_CPP_FILES src/svgfill.cpp src/arrange_polygons.cpp)

add_library(svgfill ${LIB_H_FILES} ${LIB_CPP_FILES})
target_link_libraries(svgfill Boost::boost svgpp::svgpp LibXml2::LibXml2 CGAL::CGAL)

add_executable(svgfill_exe src/main.cpp)
target_link_libraries(svgfill_exe svgfill)
set_property(TARGET svgfill_exe PROPERTY OUTPUT_NAME svgfill)
if(WIN32)
	# both the library and the executable now result in a file with basename svgfill,
	# on linux the library is prefixed with lib as libsvgfill.a. Windows does not
	# have this mechanism, so on windows the linker would be created an import library
	# for the executable, also named svgfill.lib. This naming conflict results in:
	#     LINK : fatal error LNK1149: output filename matches input filename
	# This flag tells the linker not to generate an import library and therefore no
	# conflict occurs.
    target_link_options(svgfill_exe PRIVATE "/NOIMPLIB")
endif()

set_target_properties(svgfill PROPERTIES PUBLIC_HEADER "src/svgfill.h;src/progress.h")
install(TARGETS svgfill)

install(TARGETS svgfill_exe	 DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )