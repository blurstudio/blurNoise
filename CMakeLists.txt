cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)

project(blurNoise)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
find_package(Maya REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/../output/maya/${MAYA_VERSION}/blur)


set(SOURCE_FILES
    "src/blurNoise.h"
    "src/blurNoise.cpp"
    "src/pluginMain.cpp"
    "src/pluginMain.cpp"
    "src/OpenSimplex2F.h"
    "src/OpenSimplex2F.cpp"
)

include_directories(${MAYA_INCLUDE_DIR})
link_directories(${MAYA_LIBRARY_DIR})
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${MAYA_LIBRARIES} ${OPENGL_LIBRARY} ${OpenMP_CXX_LIBRARIES})

MAYA_PLUGIN(${PROJECT_NAME})
install(TARGETS ${PROJECT_NAME} ${MAYA_TARGET_TYPE} DESTINATION plug-ins)

