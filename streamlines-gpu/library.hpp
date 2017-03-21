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

#pragma once

/* ************************************************************************ */

#ifndef OPENCL_STATIC

/* ************************************************************************ */

// OpenCL
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include "CL/opencl.h"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines_gpu {

/* ************************************************************************ */

// Function pointers
#define CL_FUNCTION(name) extern decltype(name)* name ## _fn;
#include "functions.hpp"

/* ************************************************************************ */

/**
 * @brief      Initialize OpenCL library.
 */
void library_init();

/* ************************************************************************ */

/**
 * @brief      Terminate OpenCL library.
 */
void library_free();

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */

#endif

/* ************************************************************************ */
