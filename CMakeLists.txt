cmake_minimum_required(VERSION 2.8)
project(envmap)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_STANDARD 11)

include_directories( "src/")
include_directories(SYSTEM "deps/")

add_executable(envmap src/envmap.cpp src/envUtils.cpp
  src/computeSphericalHarmonicsFromCubemap.cpp
  src/prefilterCubemapGGX.cpp
  src/utils.cpp
  src/log.cpp
  src/Image.cpp src/Cubemap.cpp )

target_link_libraries(envmap pthread)
target_compile_options(envmap PRIVATE -Wall -Wextra -pedantic -Werror)

install(TARGETS envmap
RUNTIME DESTINATION bin
)

