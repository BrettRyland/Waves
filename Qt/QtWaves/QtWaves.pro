TEMPLATE = app
TARGET = QtWaves
QT += core widgets gui
CONFIG += release console optimize_full
DEFINES += QT_DLL QT_WIDGETS_LIB CL_TARGET_OPENCL_VERSION=220
QMAKE_CXXFLAGS += -Wno-unused-variable -Wno-unused-parameter
INCLUDEPATH += ./GeneratedFiles \
    . \
    ./GeneratedFiles/Release \
    include \
    ../../include
LIBS += -lOpenCL
DEPENDPATH += .
MOC_DIR += ./GeneratedFiles/release
OBJECTS_DIR += release
UI_DIR += ./GeneratedFiles
RCC_DIR += ./GeneratedFiles
include(QtWaves.pri)
