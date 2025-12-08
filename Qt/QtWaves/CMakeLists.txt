cmake_minimum_required(VERSION 3.16)
project(QtWaves VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt6 REQUIRED COMPONENTS Gui Widgets OpenGL OpenGLWidgets)
set(CMAKE_AUTOUIC_SEARCH_PATHS forms)
# qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

add_compile_options(-fPIC -Wall -Werror -Wextra -Wno-unused-parameter -Wno-unused-variable -pedantic)

file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB INCLUDES CONFIGURE_DEPENDS include/*.h)
qt_wrap_cpp(SOURCES ${INCLUDES} TARGET QtWaves)
qt_add_executable(QtWaves ${SOURCES} ${INCLUDES})

target_include_directories(QtWaves PRIVATE include)

target_compile_definitions(QtWaves PRIVATE
    CL_TARGET_OPENCL_VERSION=300
)

target_link_libraries(QtWaves PRIVATE
    OpenCL
    Qt::Core
    Qt::Gui
    Qt::Widgets
    Qt::OpenGL
    Qt::OpenGLWidgets
)

qt_add_resources(QtWaves "qtwaves"
    PREFIX
        "/QtWaves"
    FILES
        "shaders/surface.frag.glsl"
        "shaders/surface.vert.glsl"
)

