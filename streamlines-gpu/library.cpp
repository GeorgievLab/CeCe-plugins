/* ************************************************************************ */
/* Georgiev Lab (c) 2015-2017                                               */
/* ************************************************************************ */
/* Department of Cybernetics                                                */
/* Faculty of Applied Sciences                                              */
/* University of West Bohemia in Pilsen                                     */
/* ************************************************************************ */
/*                                                                          */
/* This file is part of CeCe.                                               */
/*                                                                          */
/* CeCe is free software: you can redistribute it and/or modify             */
/* it under the terms of the GNU General Public License as published by     */
/* the Free Software Foundation, either version 3 of the License, or        */
/* (at your option) any later version.                                      */
/*                                                                          */
/* CeCe is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU General Public License for more details.                             */
/*                                                                          */
/* You should have received a copy of the GNU General Public License        */
/* along with CeCe.  If not, see <http://www.gnu.org/licenses/>.            */
/*                                                                          */
/* ************************************************************************ */

#ifndef OPENCL_STATIC

/* ************************************************************************ */

// Declaration
#include "library.hpp"

// OpenCL
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include "CL/opencl.h"

#if _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// CeCe
#include "cece/core/String.hpp"
#include "cece/core/Exception.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines_gpu {

/* ************************************************************************ */

// Function pointers
#define CL_FUNCTION(name) decltype(name)* name ## _fn;
#include "functions.hpp"

/* ************************************************************************ */

#ifndef OPENCL_STATIC
#if _WIN32
static HMODULE g_library = nullptr;
#else
static void* g_library = nullptr;
#endif
#endif

/* ************************************************************************ */

void library_init()
{
#if __linux__
    g_library = dlopen("libOpenCL.so", RTLD_LAZY);

    if (g_library == nullptr)
        g_library = dlopen("libOpenCL.so.1", RTLD_LAZY);

#elif __APPLE__ && __MACH__
    g_library = dlopen("libOpenCL.dylib", RTLD_LAZY);
#elif _WIN32
    // Disable error reporting box (on Windows 7)
    DWORD oldMode;
    SetThreadErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX, &oldMode);
    g_library = LoadLibraryW(L"OpenCL.dll");
    SetThreadErrorMode(oldMode, NULL);
#endif

#if _WIN32
    if (g_library == nullptr)
    {
        // Get error message
        const auto err = ::GetLastError();
        if (err == 0)
            return;

        LPSTR buffer = nullptr;
        auto size = FormatMessageA(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL, err, MAKELANGID(LANG_ENGLISH, SUBLANG_ENGLISH_US), (LPSTR)&buffer, 0, NULL
        );

        // Remove trailing new lines.
        while (size > 0 && (buffer[size - 1] == '\r' || buffer[size - 1] == '\n'))
        {
            --size;
        }

        String message(buffer, size);

        // Free the buffer
        LocalFree(buffer);

        throw RuntimeException(message);
    }
#else
    if (g_library == nullptr)
        throw RuntimeException(dlerror());
#endif

    // Assign
#if _WIN32
# define CL_FUNCTION(name) name ## _fn = reinterpret_cast<decltype(name)*>(GetProcAddress(g_library, # name));
#else
# define CL_FUNCTION(name) name ## _fn = reinterpret_cast<decltype(name)*>(dlsym(g_library, # name));
#endif
#include "functions.hpp"


    // Check
#define CL_FUNCTION(name) if (name ## _fn == nullptr) throw RuntimeException("Can't obtain " # name " function");
#include "functions.hpp"

}

/* ************************************************************************ */

void library_free()
{
#if _WIN32
    FreeLibrary(g_library);
#else
    if (g_library)
        dlclose(g_library);
#endif
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */

#endif

/* ************************************************************************ */
