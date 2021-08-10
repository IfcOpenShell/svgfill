cmake_minimum_required (VERSION 2.8.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (svgfill)

set(CGAL_LIBRARY_NAMES libCGAL_Core libCGAL_ImageIO libCGAL)

if("${CGAL_INCLUDE_DIR}" STREQUAL "")
	set(CGAL_INCLUDE_DIR "/usr/include/" CACHE FILEPATH "CGAL header files")
	message(STATUS "Looking for CGAL include files in: ${CGAL_INCLUDE_DIR}")
	message(STATUS "Use CGAL_INCLUDE_DIR to specify another directory")
else()
	set(CGAL_INCLUDE_DIR ${CGAL_INCLUDE_DIR} CACHE FILEPATH "CGAL header files")
	message(STATUS "Looking for CGAL include files in: ${CGAL_INCLUDE_DIR}")
endif()

if(NOT "${CGAL_LIBRARY_DIR}" STREQUAL "")
	set(CGAL_LIBRARY_DIR ${CGAL_LIBRARY_DIR} CACHE FILEPATH "CGAL library files")
	message(STATUS "Looking for CGAL library files in: ${CGAL_LIBRARY_DIR}")
endif()

find_library(libCGAL NAMES CGAL PATHS ${CGAL_LIBRARY_DIR} NO_DEFAULT_PATH)

if(libCGAL)
	message(STATUS "CGAL library files found")
	foreach(lib ${CGAL_LIBRARY_NAMES})
		string(REPLACE libCGAL "${lib}" lib_path "${libCGAL}")
		list(APPEND CGAL_LIBRARIES "${lib_path}")
	endforeach()
else()
	if(NOT "${CGAL_LIBRARY_DIR}" STREQUAL "")
		file(GLOB CGAL_LIBRARIES ${CGAL_LIBRARY_DIR}/CGAL*.lib)
		list(LENGTH CGAL_LIBRARY_NAMES num_cgal_library_names)
		list(LENGTH CGAL_LIBRARIES num_cgal_libraries)
		link_directories("${CGAL_LIBRARY_DIR}")
		if(NOT "${num_cgal_library_names}" STREQUAL "${num_cgal_libraries}")
			message(FATAL_ERROR "Unable to find CGAL library files, aborting")
		endif()
		message(STATUS "CGAL library files found")
	endif()
endif()

find_library(libGMP NAMES libgmp.a gmp mpir PATHS ${GMP_LIBRARY_DIR} NO_DEFAULT_PATH)
find_library(libMPFR NAMES libmpfr.a mpfr PATHS ${MPFR_LIBRARY_DIR} NO_DEFAULT_PATH)
if(NOT libGMP)
        message(FATAL_ERROR "Unable to find GMP library files, aborting")
endif()
if(NOT libMPFR)
        message(FATAL_ERROR "Unable to find MPFR library files, aborting")
endif()

list(APPEND CGAL_LIBRARIES "${libMPFR}")
list(APPEND CGAL_LIBRARIES "${libGMP}")

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

find_package(Boost)
message(STATUS "Boost include files found in ${Boost_INCLUDE_DIRS}")

include_directories(${Boost_INCLUDE_DIRS} ${LIBXML2_INCLUDE_DIR}
	${CGAL_INCLUDE_DIR} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/svgpp/include
)

file(GLOB LIB_H_FILES src/*.h)
file(GLOB LIB_CPP_FILES src/svgfill.cpp)
set(LIB_SRC_FILES ${LIB_H_FILES} ${LIB_CPP_FILES})
add_library(libsvgfill ${LIB_SRC_FILES})
target_link_libraries(libsvgfill ${Boost_LIBRARIES} ${BCRYPT_LIBRARIES} ${LIBXML2_LIBRARIES} ${CGAL_LIBRARIES})

add_executable(svgfill src/main.cpp)
target_link_libraries(svgfill libsvgfill)

install(TARGETS svgfill DESTINATION bin)
install(TARGETS libsvgfill DESTINATION lib)
install(FILES ${LIB_H_FILES} DESTINATION include)
