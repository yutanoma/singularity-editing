cmake_minimum_required(VERSION 3.16)
project(singularity-editing)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   ON)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
# option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
# option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
# option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)
# option(LIBIGL_WITH_XML               "Use XML"                      ON)

# directional
include(directional)

include(libigl)

include(geometry-central)

set(DATA_PATH ${CMAKE_CURRENT_LIST_DIR}/data CACHE PATH "location of data path")

# Add your project files
FILE(GLOB_RECURSE SRC ./src/*.cc ./src/*.cpp)
# aux_source_directory(utils UTL_SRC)
# aux_source_directory(models MODEL_SRC)
# aux_source_directory(controllers CTRL_SRC)
# aux_source_directory(views VIEW_SRC)
# aux_source_directory(services SRV_SRC)

add_executable(${PROJECT_NAME} main.cc)

add_library(env INTERFACE)
target_compile_definitions(env INTERFACE "-DDATA_PATH=\"${DATA_PATH}\"")

target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui igl::triangle directional saddlepoint env geometry-central)
target_sources(${PROJECT_NAME} PUBLIC ${SRC})
# target_sources(${PROJECT_NAME} PUBLIC ${UTL_SRC} ${MODEL_SRC} ${CTRL_SRC} ${VIEW_SRC} ${SRV_SRC})
# target_include_directories(${PROJECT_NAME} PUBLIC views)
