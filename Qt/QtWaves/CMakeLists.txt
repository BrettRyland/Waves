cmake_minimum_required(VERSION 3.16)
project(QtWaves VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(Qt6_DIR "C:/Qt/6.12.0/msvc2022_64/lib/cmake/Qt6")
find_package(Qt6 REQUIRED COMPONENTS Gui Widgets OpenGL OpenGLWidgets)
set(CMAKE_AUTOUIC_SEARCH_PATHS forms)
# qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS compute)
find_package(OpenCL REQUIRED)
find_package(OpenGL REQUIRED)

if (MSVC)
    add_compile_options(/W4 /WX /wd4100 /wd4101 /wd4189 /wd4458 /permissive- /EHsc)
else()
    add_compile_options(-fPIC -Wall -Werror -Wextra -Wno-unused-parameter -Wno-unused-variable -pedantic)
endif()

file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB INCLUDES CONFIGURE_DEPENDS include/*.h)
qt_wrap_cpp(SOURCES ${INCLUDES} TARGET QtWaves)
qt_add_executable(QtWaves ${SOURCES} ${INCLUDES})

target_include_directories(QtWaves PRIVATE
    include
    ${Boost_INCLUDE_DIRS}
    ${OpenCL_INCLUDE_DIRS}
)

target_compile_definitions(QtWaves PRIVATE
    CL_TARGET_OPENCL_VERSION=300
)

target_link_libraries(QtWaves PRIVATE
    OpenCL::OpenCL
    OpenGL::GL
    Boost::compute
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

if (WIN32)
	set(SHADER_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/shaders")
	set(KERNEL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kernels")
	set(SHADER_DEST_DIR "$<TARGET_FILE_DIR:QtWaves>/shaders")
	set(KERNEL_DEST_DIR "$<TARGET_FILE_DIR:QtWaves>/kernels")
    add_custom_command(
        TARGET QtWaves
        POST_BUILD

        COMMENT "Copying OpenCL kernels and shaders to build output directory: from ${SHADER_SOURCE_DIR} to ${SHADER_DEST_DIR}, from ${KERNEL_SOURCE_DIR} to ${KERNEL_DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E make_directory ${SHADER_DEST_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${SHADER_SOURCE_DIR} ${SHADER_DEST_DIR}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${KERNEL_DEST_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${KERNEL_SOURCE_DIR} ${KERNEL_DEST_DIR}

        COMMENT "Running windeployqt to include the appropriate libraries."
        COMMAND "C:/Qt/6.12.0/msvc2022_64/bin/windeployqt.exe" "$<TARGET_FILE:QtWaves>"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endif()