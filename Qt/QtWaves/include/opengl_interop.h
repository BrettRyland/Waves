#ifdef _WIN32
#define NOMINMAX // Otherwise windows.h breaks boost.compute min/max definitions
#include <windows.h>
#include <wingdi.h> // Contains wglGetCurrentContext
#endif
#include <boost/compute/interop/opengl.hpp>