include(FetchContent)

## include saddlepoint

FetchContent_Declare(
    saddlepoint
    GIT_REPOSITORY https://github.com/avaxman/SaddlePoint.git
    GIT_TAG 4f95623
)

FetchContent_GetProperties(saddlepoint)
if(NOT saddlepoint_POPULATED)
    FetchContent_Populate(saddlepoint)
endif()

add_library(saddlepoint INTERFACE)
target_include_directories(saddlepoint INTERFACE ${saddlepoint_SOURCE_DIR}/include)

## include directional

FetchContent_Declare(
    directional
    GIT_REPOSITORY https://github.com/avaxman/Directional.git
    GIT_TAG 46d8d3c
)

FetchContent_GetProperties(directional)
if(NOT directional_POPULATED)
    FetchContent_Populate(directional)
endif()

add_library(directional INTERFACE)
target_include_directories(directional INTERFACE ${directional_SOURCE_DIR}/include)
