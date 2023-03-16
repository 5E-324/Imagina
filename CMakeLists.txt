cmake_minimum_required (VERSION 3.8)

project ("Imagina")

option(ENABLE_CUSTOM_FORMULA "Enable custom formula" OFF)
option(MPIR_USE_DLL "Use dll version of MPIR" ON)
set(INCLUDE_PATH "include" CACHE STRING "Include path")
set(LIBRARY_PATH "lib" CACHE STRING "Library path")

FILE(GLOB SourceFiles "*.cpp")

if(ENABLE_CUSTOM_FORMULA)
	FILE(GLOB SourceFilesEncoder "Encoder/*.cpp")
else()
	set(SourceFilesEncoder "")
endif()

add_executable (Imagina WIN32 ${SourceFiles} ${SourceFilesEncoder} "Resource.rc")
target_precompile_headers (Imagina PRIVATE "Includes.h")

set_source_files_properties("Resource.rc" PROPERTIES LANGUAGE RC)

target_compile_definitions(Imagina PRIVATE UNICODE)

target_include_directories(Imagina PRIVATE SYSTEM ${INCLUDE_PATH})
target_link_directories(Imagina PRIVATE ${LIBRARY_PATH})
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_include_directories(Imagina PRIVATE SYSTEM "${INCLUDE_PATH}/Debug")
    target_link_directories(Imagina PRIVATE "${LIBRARY_PATH}/Debug")
else()
	target_include_directories(Imagina PRIVATE SYSTEM "${INCLUDE_PATH}/Release")
    target_link_directories(Imagina PRIVATE "${LIBRARY_PATH}/Release")
endif()

target_link_libraries(Imagina "opengl32")
target_link_libraries(Imagina "comctl32")
target_link_libraries(Imagina "libpng16")
target_link_libraries(Imagina "zlibwapi")
target_link_libraries(Imagina "dwmapi")
target_link_libraries(Imagina "mpir")

if(MPIR_USE_DLL)
	target_compile_definitions(Imagina PRIVATE MSC_USE_DLL)
else()
	target_link_libraries(Imagina "mpirxx")
endif()

if(ENABLE_CUSTOM_FORMULA)
	target_link_libraries(Imagina "mpfr")
	target_link_libraries(Imagina "mpc")
	target_link_libraries(Imagina "symengine")
	
	target_compile_definitions(Imagina PRIVATE ENABLE_CUSTOM_FORMULA)
endif()

set_property(TARGET Imagina PROPERTY CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno")

set_source_files_properties(HInfLAEvaluator.cpp PROPERTIES COMPILE_FLAGS "-fassociative-math -ffp-contract=fast" SKIP_PRECOMPILE_HEADERS ON)

set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -Xlinker /MANIFEST:NO")
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -Xlinker /ignore:4099")