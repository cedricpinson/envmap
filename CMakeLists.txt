cmake_minimum_required(VERSION 2.8)
project(envmap)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

include_directories( "src/")
include_directories(SYSTEM "deps/")

set(SOURCES
  src/envmap.cpp src/envUtils.cpp
  src/computeSphericalHarmonicsFromCubemap.cpp
  src/environmentConversion.cpp
  src/prefilterCubemapGGX.cpp
  src/package.cpp
  src/utils.cpp
  src/io.cpp
  src/log.cpp
  src/encode.cpp
  src/Image.cpp
  src/Cubemap.cpp
  src/SummedAreaTable.cpp
  src/lightExtraction.cpp
)

set(CLANG_FORMAT_EXCLUDE_PATTERNS "build/" "debug/" "deps/")
include("${CMAKE_MODULE_PATH}/ClangFormat.cmake")

add_executable(envmap ${SOURCES})

target_link_libraries(envmap pthread)
target_compile_options(envmap PRIVATE -Wall -Wextra -pedantic -Werror -Wdouble-promotion -Wmissing-variable-declarations -Wfloat-conversion -Wpadded -g)

install(TARGETS envmap
  RUNTIME DESTINATION bin
)

add_executable(sampling
  tests/sampling.cpp
  src/Cubemap.cpp
  src/Image.cpp
  )

target_link_libraries(sampling)
target_compile_options(sampling PRIVATE -Wall -Wextra -pedantic -Werror -Wdouble-promotion -Wmissing-variable-declarations -Wfloat-conversion -Wpadded -g)

add_custom_target(check)

add_custom_command( TARGET check
  POST_BUILD
  COMMAND ./sampling
)
