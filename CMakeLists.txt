cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(thorax_truetype VERSION 0.0.1)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(thorax_truetype "source/thorax_truetype.cpp"
	"source/renderer/FontRenderInfo.cpp"
	"source/renderer/ShapeBuilder.cpp"
	"source/renderer/bvh.cpp"
	"source/renderer/drawcontext.cpp" 
	"source/grid_builder.cpp"
	"source/thorax_bridge.cpp"
	)

if(MSVC)
    target_compile_options(thorax_truetype PRIVATE /arch:AVX2)
endif()

target_include_directories(thorax_truetype PUBLIC "include")
