cmake_minimum_required(VERSION 2.8)
project(envmap)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_STANDARD 11)

find_package(TBB)

include_directories( "deps/" "src/")

add_executable(envmap src/envmap.cpp src/envUtils.cpp
  src/computeSphericalHarmonicsFromCubemap.cpp
  src/prefilterCubemapGGX.cpp
  src/utils.cpp
  src/Image.cpp src/Cubemap.cpp )

target_link_libraries(envmap pthread)

install(TARGETS envmap
RUNTIME DESTINATION bin
)

