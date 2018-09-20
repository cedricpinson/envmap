cmake_minimum_required(VERSION 2.8)
project(envmap)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_STANDARD 11)

include_directories( "src/")
include_directories(SYSTEM "deps/")

add_executable(envmap src/envmap.cpp src/envUtils.cpp
  src/computeSphericalHarmonicsFromCubemap.cpp
  src/environmentConversion.cpp
  src/prefilterCubemapGGX.cpp
  src/package.cpp
  src/utils.cpp
  src/io.cpp
  src/log.cpp
  src/encode.cpp
  src/Image.cpp src/Cubemap.cpp )

target_link_libraries(envmap pthread)
target_compile_options(envmap PRIVATE -Wall -Wextra -pedantic -Werror -Wdouble-promotion -Wmissing-variable-declarations -Wfloat-conversion -Wpadded -g)

install(TARGETS envmap
RUNTIME DESTINATION bin
)

