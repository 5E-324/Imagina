cmake_minimum_required (VERSION 3.8)

option(EnableCustomFormula "EnableCustomFormula" OFF)
message(${EnableCustomFormula})
FILE(GLOB SourceFiles "*.cpp")

if(EnableCustomFormula)
	FILE(GLOB SourceFilesEncoder "Encoder/*.cpp")
else()
	set(SourceFilesEncoder "")
endif()

add_executable (Imagina WIN32 ${SourceFiles} ${SourceFilesEncoder} "Resource.rc")
target_precompile_headers (Imagina PRIVATE "Includes.h")

set_source_files_properties("Resource.rc" PROPERTIES LANGUAGE RC)

target_compile_definitions(Imagina PRIVATE UNICODE)

target_include_directories(Imagina PRIVATE ".")
include_directories(SYSTEM "include")
target_include_directories(Imagina PRIVATE "C:/Program Files (x86)/Windows Kits/10/Include/10.0.19041.0/um")
target_link_directories(Imagina PRIVATE "lib/x64")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    target_link_directories(Imagina PRIVATE "lib/x64/Debug")
	include_directories(SYSTEM "include/Debug")
else()
    target_link_directories(Imagina PRIVATE "lib/x64/Release")
	include_directories(SYSTEM "include/Release")
endif()

target_link_libraries(Imagina "opengl32")
target_link_libraries(Imagina "mpir")
target_link_libraries(Imagina "comctl32")
target_link_libraries(Imagina "libpng16")
target_link_libraries(Imagina "zlibwapi")
target_link_libraries(Imagina "dwmapi")

if(EnableCustomFormula)
	target_link_libraries(Imagina "mpfr")
	target_link_libraries(Imagina "mpc")
	target_link_libraries(Imagina "symengine")
	
	target_compile_definitions(Imagina PRIVATE ENABLE_CUSTOM_FORMULA)
endif()

set_property(TARGET Imagina PROPERTY CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -DMSC_USE_DLL")

set_source_files_properties(HInfLAEvaluator.cpp PROPERTIES COMPILE_FLAGS "-fassociative-math -ffp-contract=fast" SKIP_PRECOMPILE_HEADERS ON)

set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -Xlinker /MANIFEST:NO")
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -Xlinker /ignore:4099")