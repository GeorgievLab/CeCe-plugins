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

#ifndef CL_FUNCTION
#define CL_FUNCTION(name)
#endif

/* ************************************************************************ */

CL_FUNCTION(clGetPlatformInfo)
CL_FUNCTION(clGetPlatformIDs)
CL_FUNCTION(clGetDeviceIDs)
CL_FUNCTION(clCreateContext)
CL_FUNCTION(clReleaseContext)
CL_FUNCTION(clCreateCommandQueue)
CL_FUNCTION(clReleaseCommandQueue)
CL_FUNCTION(clCreateProgramWithSource)
CL_FUNCTION(clReleaseProgram)
CL_FUNCTION(clBuildProgram)
CL_FUNCTION(clGetProgramBuildInfo)
CL_FUNCTION(clCreateKernel)
CL_FUNCTION(clReleaseKernel)
CL_FUNCTION(clCreateBuffer)
CL_FUNCTION(clReleaseMemObject)
CL_FUNCTION(clFinish)
CL_FUNCTION(clEnqueueReadBuffer)
CL_FUNCTION(clSetKernelArg)
CL_FUNCTION(clEnqueueNDRangeKernel)
CL_FUNCTION(clEnqueueBarrierWithWaitList)
CL_FUNCTION(clEnqueueWriteBuffer)
CL_FUNCTION(clReleaseEvent)

/* ************************************************************************ */

#undef CL_FUNCTION

/* ************************************************************************ */
