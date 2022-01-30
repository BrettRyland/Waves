TEMPLATE = app
TARGET = QtWaves
QT += core widgets gui
CONFIG += debug console
DEFINES += QT_DLL QT_WIDGETS_LIB
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
