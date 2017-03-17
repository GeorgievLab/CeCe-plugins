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

// Declaration
#include "Lattice.hpp"

// CeCe
#include "cece/core/DynamicArray.hpp"
#include "cece/core/VectorRange.hpp"
#include "cece/core/Exception.hpp"
#include "cece/core/String.hpp"
#include "cece/core/Log.hpp"

// Plugin
#include "../Descriptor.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace gpu {

/* ************************************************************************ */

#define CL_CHECK_RET(ret) if (ret != CL_SUCCESS) throw RuntimeException(clMessage(ret))
#define CL_CHECK(cmd) do { const cl_int ret = cmd; CL_CHECK_RET(ret); } while (false)

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

#include "program.hpp"

/* ************************************************************************ */

/**
 * @brief      Returns string representation of OpenCL error code.
 *
 * @param[in]  err   The error code.
 *
 * @return     String representation.
 */
const char* clMessage(cl_int err) noexcept
{
#define CASE(code, msg) case code: return msg;

    switch (err)
    {
    default: return "Unknown";
    CASE(CL_SUCCESS,                                    "Success")
    CASE(CL_DEVICE_NOT_FOUND,                           "Device not found")
    CASE(CL_DEVICE_NOT_AVAILABLE,                       "Device not available")
    CASE(CL_COMPILER_NOT_AVAILABLE,                     "Compiler notavailable")
    CASE(CL_MEM_OBJECT_ALLOCATION_FAILURE,              "Mem object allocation failure")
    CASE(CL_OUT_OF_RESOURCES,                           "Out of resources")
    CASE(CL_OUT_OF_HOST_MEMORY,                         "Out of host_memory")
    CASE(CL_PROFILING_INFO_NOT_AVAILABLE,               "Profiling info not available")
    CASE(CL_MEM_COPY_OVERLAP,                           "Mem copy overlap")
    CASE(CL_IMAGE_FORMAT_MISMATCH,                      "Image format mismatch")
    CASE(CL_IMAGE_FORMAT_NOT_SUPPORTED,                 "Image format not supported")
    CASE(CL_BUILD_PROGRAM_FAILURE,                      "Build program failure")
    CASE(CL_MAP_FAILURE,                                "Map failure")
    CASE(CL_MISALIGNED_SUB_BUFFER_OFFSET,               "Misaligned sub buffer offset")
    CASE(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST,  "Exec status error for events in wait list")
    CASE(CL_COMPILE_PROGRAM_FAILURE,                    "Compile program failure")
    CASE(CL_LINKER_NOT_AVAILABLE,                       "Linker not available")
    CASE(CL_LINK_PROGRAM_FAILURE,                       "Link program failure")
    CASE(CL_DEVICE_PARTITION_FAILED,                    "Device partition failed")
    CASE(CL_KERNEL_ARG_INFO_NOT_AVAILABLE,              "Kernel arg info not available")
    CASE(CL_INVALID_VALUE,                              "Invalid value")
    CASE(CL_INVALID_DEVICE_TYPE,                        "Invalid device type")
    CASE(CL_INVALID_PLATFORM,                           "Invalid platform")
    CASE(CL_INVALID_DEVICE,                             "Invalid device")
    CASE(CL_INVALID_CONTEXT,                            "Invalid context")
    CASE(CL_INVALID_QUEUE_PROPERTIES,                   "Invalid queue properties")
    CASE(CL_INVALID_COMMAND_QUEUE,                      "Invalid command queue")
    CASE(CL_INVALID_HOST_PTR,                           "Invalid host ptr")
    CASE(CL_INVALID_MEM_OBJECT,                         "Invalid mem object")
    CASE(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR,            "Invalid image format descriptor")
    CASE(CL_INVALID_IMAGE_SIZE,                         "Invalid image_size")
    CASE(CL_INVALID_SAMPLER,                            "Invalid sampler")
    CASE(CL_INVALID_BINARY,                             "Invalid binary")
    CASE(CL_INVALID_BUILD_OPTIONS,                      "Invalid build options")
    CASE(CL_INVALID_PROGRAM,                            "Invalid program")
    CASE(CL_INVALID_PROGRAM_EXECUTABLE,                 "Invalid program executable")
    CASE(CL_INVALID_KERNEL_NAME,                        "Invalid kernel name")
    CASE(CL_INVALID_KERNEL_DEFINITION,                  "Invalid kernel definition")
    CASE(CL_INVALID_KERNEL,                             "Invalid kernel")
    CASE(CL_INVALID_ARG_INDEX,                          "Invalid arg index")
    CASE(CL_INVALID_ARG_VALUE,                          "Invalid arg value")
    CASE(CL_INVALID_ARG_SIZE,                           "Invalid arg size")
    CASE(CL_INVALID_KERNEL_ARGS,                        "Invalid kernel args")
    CASE(CL_INVALID_WORK_DIMENSION,                     "Invalid work dimension")
    CASE(CL_INVALID_WORK_GROUP_SIZE,                    "Invalid work group size")
    CASE(CL_INVALID_WORK_ITEM_SIZE,                     "Invalid work item size")
    CASE(CL_INVALID_GLOBAL_OFFSET,                      "Invalid global offset")
    CASE(CL_INVALID_EVENT_WAIT_LIST,                    "Invalid event wait_list")
    CASE(CL_INVALID_EVENT,                              "Invalid event")
    CASE(CL_INVALID_OPERATION,                          "Invalid operation")
    CASE(CL_INVALID_GL_OBJECT,                          "Invalid GL object")
    CASE(CL_INVALID_BUFFER_SIZE,                        "Invalid buffer size")
    CASE(CL_INVALID_MIP_LEVEL,                          "Invalid mip level")
    CASE(CL_INVALID_GLOBAL_WORK_SIZE,                   "Invalid global work size")
    CASE(CL_INVALID_PROPERTY,                           "Invalid property")
    CASE(CL_INVALID_IMAGE_DESCRIPTOR,                   "Invalid image descriptor")
    CASE(CL_INVALID_COMPILER_OPTIONS,                   "Invalid compiler options")
    CASE(CL_INVALID_LINKER_OPTIONS,                     "Invalid linker options")
    CASE(CL_INVALID_DEVICE_PARTITION_COUNT,             "Invalid device partition count")
#if CL_VERSION_2_0
    CASE(CL_INVALID_PIPE_SIZE,                          "Invalid pipe size")
    CASE(CL_INVALID_DEVICE_QUEUE,                       "Invalid device queue")
#endif
    }
}

/* ************************************************************************ */

/**
 * @brief      Get platform info.
 *
 * @param[in]  platform  The platform.
 *
 * @return     The platform information.
 */
std::string getPlatformInfo(cl_platform_id platform, cl_platform_info info)
{
    size_t size;
    CL_CHECK(clGetPlatformInfo(platform, info, 0, nullptr, &size));

    std::vector<char> buffer(size);

    CL_CHECK(clGetPlatformInfo(platform, info, size, buffer.data(), nullptr));

    return String(buffer.data(), size);
}

/* ************************************************************************ */

/**
 * @brief      Choose best OpenCL platform.
 *
 * @return     Platform ID.
 */
cl_platform_id choosePlatform()
{
    // Get first platform
    cl_platform_id platform;
    CL_CHECK(clGetPlatformIDs(1, &platform, nullptr));

    // Print platform info
    Log::info("[streamlines] OpenCL platform name: ", getPlatformInfo(platform, CL_PLATFORM_NAME));
    Log::info("[streamlines] OpenCL platform vendor: ", getPlatformInfo(platform, CL_PLATFORM_VENDOR));
    Log::info("[streamlines] OpenCL platform version: ", getPlatformInfo(platform, CL_PLATFORM_VERSION));
    Log::info("[streamlines] OpenCL platform profile: ", getPlatformInfo(platform, CL_PLATFORM_PROFILE));

    return platform;
}

/* ************************************************************************ */

/**
 * @brief      Choose best OpenCL device.
 *
 * @param[in]  platform  The platform.
 *
 * @return     Device.
 */
cl_device_id chooseDevice(cl_platform_id platform)
{
    // Get default device
    cl_device_id device;
    CL_CHECK(clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr));

    return device;
}

/* ************************************************************************ */

/**
 * @brief      Error callback.
 *
 * @param[in]  error         The error.
 * @param[in]  private_info  The private information.
 * @param[in]  cb            The private information size.
 * @param      user_data     The user data.
 */
void error_callback(const char* error, const void* private_info, size_t cb, void* user_data)
{
    Log::error("[streamlines] OpenCL error: ", error);
}

/* ************************************************************************ */

/**
 * @brief      Calculates the offset in 2D array.
 *
 * @param[in]  size  The size.
 * @param[in]  x     The X coordinate.
 * @param[in]  y     The Y coordinate.
 *
 * @return     The array index/offset.
 */
size_t calcOffset(Lattice::SizeType size, size_t x, size_t y) noexcept
{
    return y * size.getWidth() + x;
}

/* ************************************************************************ */

/**
 * @brief      Calculates the offset in "3D" distribution function array.
 *
 * @param[in]  size  The size.
 * @param[in]  x     The X coordinate.
 * @param[in]  y     The Y coordinate.
 * @param[in]  iPop  The pop
 *
 * @return     The array index/offset.
 */
size_t calcDfOffset(Lattice::SizeType size, size_t x, size_t y, int iPop) noexcept
{
    return size.getWidth() * size.getHeight() * iPop + calcOffset(size, x, y);
}

/* ************************************************************************ */

/**
 * @brief      Convert coordinate to BC position.
 *
 * @return     Number of BC position.
 */
cl_int coordToBcPosition(Lattice::SizeType size, Lattice::CoordinateType coord)
{
    if (coord.getX() == 0)
        return static_cast<cl_int>(Lattice::BoundaryPosition::Left);
    else if (coord.getX() == size.getX() - 1)
        return static_cast<cl_int>(Lattice::BoundaryPosition::Right);
    else if (coord.getY() == 0)
        return static_cast<cl_int>(Lattice::BoundaryPosition::Top);
    else if (coord.getY() == size.getY() - 1)
        return static_cast<cl_int>(Lattice::BoundaryPosition::Bottom);

    throw RuntimeException("Unknown boundary position");
}

/* ************************************************************************ */

}

/* ************************************************************************ */

Lattice::Lattice(SizeType size, RealType omega)
    : streamlines::Lattice(size, omega)
    , m_platform(choosePlatform())
    , m_device(chooseDevice(m_platform))
{
    // Create OpenCL context
    {
        cl_int ret;
        m_context = clCreateContext(nullptr, 1, &m_device, error_callback, nullptr, &ret);
        CL_CHECK_RET(ret);
    }

    // Create command queue
    {
        cl_int ret;
        m_commandQueue = clCreateCommandQueue(m_context, m_device, 0, &ret);
        CL_CHECK_RET(ret);
    }

    // Create program
    {
        cl_int ret;
        const size_t length = strlen(PROGRAM);
        const char* source = PROGRAM;

        // Create program from source
        m_program = clCreateProgramWithSource(m_context, 1, &source, &length, &ret);
        CL_CHECK_RET(ret);
    }

    // Build program
    {
        const cl_int ret = clBuildProgram(m_program, 1, &m_device, nullptr, nullptr, nullptr);

        if (ret != CL_SUCCESS)
        {
            size_t len;
            clGetProgramBuildInfo(m_program, m_device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &len);

            std::string log(len, '\0');
            clGetProgramBuildInfo(m_program, m_device, CL_PROGRAM_BUILD_LOG, len, &log[0], nullptr);

            throw RuntimeException(log);
        }
    }

    // Create kernels
    {
        cl_int ret;

        m_initKernel = clCreateKernel(m_program, "init", &ret);
        CL_CHECK_RET(ret);

        m_collideKernel = clCreateKernel(m_program, "collide", &ret);
        CL_CHECK_RET(ret);

        m_streamKernel = clCreateKernel(m_program, "stream", &ret);
        CL_CHECK_RET(ret);

        m_bcKernel = clCreateKernel(m_program, "bc", &ret);
        CL_CHECK_RET(ret);

        m_syncKernel = clCreateKernel(m_program, "sync", &ret);
        CL_CHECK_RET(ret);
    }

    // Allocate memory DF
    {
        cl_int ret;
        const size_t size = Descriptor::SIZE * getSize().getWidth() * getSize().getHeight();

        // GPU buffer
        m_df = clCreateBuffer(m_context, CL_MEM_READ_WRITE, sizeof(cl_float) * size, nullptr, &ret);
        CL_CHECK_RET(ret);

        // GPU buffer
        m_df2 = clCreateBuffer(m_context, CL_MEM_READ_WRITE, sizeof(cl_float) * size, nullptr, &ret);
        CL_CHECK_RET(ret);

        // Local buffer
        m_dfLocal.buffer.resize(size);
    }

    // Allocate memory velocity
    {
        cl_int ret;
        const size_t size = getSize().getWidth() * getSize().getHeight();

        // GPU buffer
        m_velocity = clCreateBuffer(m_context, CL_MEM_READ_WRITE, sizeof(cl_float2) * size, nullptr, &ret);
        CL_CHECK_RET(ret);

        // Local buffer
        m_velocityLocal.buffer.resize(size);
    }

    // Allocate memory density
    {
        cl_int ret;
        const size_t size = getSize().getWidth() * getSize().getHeight();

        Log::debug("[streamlines] Create density buffer");

        // GPU buffer
        m_density = clCreateBuffer(m_context, CL_MEM_READ_WRITE, sizeof(cl_float) * size, nullptr, &ret);
        CL_CHECK_RET(ret);

        // Local buffer
        m_densityLocal.buffer.resize(size);
    }

    // Allocate memory Data
    {
        cl_int ret;
        const size_t size = getSize().getWidth() * getSize().getHeight();

        // Local buffer
        m_dataLocal.buffer.resize(size);

        Log::debug("[streamlines] Create data buffer");

        // GPU buffer
        m_data = clCreateBuffer(m_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Data) * size, m_dataLocal.buffer.data(), &ret);
        CL_CHECK_RET(ret);
    }
}

/* ************************************************************************ */

Lattice::~Lattice()
{
    if (m_data)
        clReleaseMemObject(m_data);

    if (m_density)
        clReleaseMemObject(m_density);

    if (m_velocity)
        clReleaseMemObject(m_velocity);

    if (m_df)
        clReleaseMemObject(m_df);

    if (m_df2)
        clReleaseMemObject(m_df2);

    if (m_syncKernel)
        clReleaseKernel(m_syncKernel);

    if (m_bcKernel)
        clReleaseKernel(m_bcKernel);

    if (m_streamKernel)
        clReleaseKernel(m_streamKernel);

    if (m_collideKernel)
        clReleaseKernel(m_collideKernel);

    if (m_initKernel)
        clReleaseKernel(m_initKernel);

    if (m_program)
        clReleaseProgram(m_program);

    if (m_commandQueue)
        clReleaseCommandQueue(m_commandQueue);

    if (m_context)
        clReleaseContext(m_context);
}

/* ************************************************************************ */

Dynamics Lattice::getDynamics(CoordinateType coord) const noexcept
{
    // Dynamics cannot be changed by GPU, so we can access directly the CPU
    // memory.

    // Calculate index
    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());
    const auto& data = m_dataLocal.buffer[i];

    return static_cast<Dynamics>(data.dynamics);
}

/* ************************************************************************ */

Lattice::VelocityType Lattice::getVelocity(CoordinateType coord) const noexcept
{
    if (m_velocityLocal.dirty)
    {
        const size_t size = getSize().getWidth() * getSize().getHeight();

        // Finish before read
        CL_CHECK(clFinish(m_commandQueue));

        // Read from GPU
        CL_CHECK(clEnqueueReadBuffer(m_commandQueue,
            m_velocity,
            CL_TRUE,
            0,
            size * sizeof(cl_float2),
            m_velocityLocal.buffer.data(),
            0, nullptr, nullptr
        ));

        m_velocityLocal.dirty = false;
    }

    // Calculate index
    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());
    const auto vel = m_velocityLocal.buffer[i];

    return {vel.s[0], vel.s[1]};
}

/* ************************************************************************ */

void Lattice::setVelocity(CoordinateType coord, VelocityType velocity) noexcept
{

}

/* ************************************************************************ */

Lattice::DensityType Lattice::getDensity(CoordinateType coord) const noexcept
{
    if (m_densityLocal.dirty)
    {
        const size_t size = getSize().getWidth() * getSize().getHeight();

        // Finish before read
        CL_CHECK(clFinish(m_commandQueue));

        // Read from GPU
        CL_CHECK(clEnqueueReadBuffer(m_commandQueue,
            m_density,
            CL_TRUE,
            0,
            size * sizeof(cl_float),
            m_densityLocal.buffer.data(),
            0, nullptr, nullptr
        ));

        m_densityLocal.dirty = false;
    }

    // Calculate index
    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());
    const auto rho = m_densityLocal.buffer[i];

    return rho;
}

/* ************************************************************************ */

void Lattice::setDensity(CoordinateType coord, DensityType density) noexcept
{

}

/* ************************************************************************ */

void Lattice::setVelocityDensity(CoordinateType coord, VelocityType velocity, DensityType density) noexcept
{

}

/* ************************************************************************ */

Lattice::DistributionsType Lattice::getDistributions(CoordinateType coord) const noexcept
{
    if (m_dfLocal.dirty)
    {
        const size_t size = Descriptor::SIZE * getSize().getWidth() * getSize().getHeight();

        // Finish before read
        CL_CHECK(clFinish(m_commandQueue));

        // Read from GPU
        CL_CHECK(clEnqueueReadBuffer(m_commandQueue,
            m_df,
            CL_TRUE,
            0,
            size * sizeof(cl_float),
            m_dfLocal.buffer.data(),
            0, nullptr, nullptr
        ));

        m_dfLocal.dirty = false;
    }

    Lattice::DistributionsType res;

    for (Descriptor::PopIndexType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        // Calculate index
        const auto i = calcDfOffset(getSize(), coord.getX(), coord.getY(), iPop);
        const auto f = m_dfLocal.buffer[i];

        // Store
        res[iPop] = f;
    }

    return res;
}

/* ************************************************************************ */

void Lattice::setDistributions(CoordinateType coord, DistributionsType distributions) noexcept
{
    CECE_ASSERT(inRange(coord));

    for (Descriptor::PopIndexType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        const auto i = calcDfOffset(getSize(), coord.getX(), coord.getY(), iPop);

        m_dfLocal.buffer[i] = distributions[iPop];
    }

    m_dfUpload = true;
}

/* ************************************************************************ */

void Lattice::setNoneDynamics(CoordinateType coord) noexcept
{
    CECE_ASSERT(inRange(coord));

    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());

    Data& data = m_dataLocal.buffer[i];
    data.dynamics = static_cast<cl_int>(Dynamics::None);

    m_dataLocal.dirty = true;
}

/* ************************************************************************ */

void Lattice::setFluidDynamics(CoordinateType coord) noexcept
{
    CECE_ASSERT(inRange(coord));

    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());

    Data& data = m_dataLocal.buffer[i];
    data.dynamics = static_cast<cl_int>(Dynamics::Fluid);

    m_dataLocal.dirty = true;
}

/* ************************************************************************ */

void Lattice::setWallDynamics(CoordinateType coord) noexcept
{
    CECE_ASSERT(inRange(coord));

    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());

    Data& data = m_dataLocal.buffer[i];
    data.dynamics = static_cast<cl_int>(Dynamics::Wall);

    m_dataLocal.dirty = true;
}

/* ************************************************************************ */

void Lattice::setInletDynamics(CoordinateType coord, VelocityType velocity) noexcept
{
    CECE_ASSERT(inRange(coord));

    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());

    Data& data = m_dataLocal.buffer[i];
    data.dynamics = static_cast<cl_int>(Dynamics::Inlet);
    data.inlet.position = coordToBcPosition(getSize(), coord);
    data.inlet.velocity = {
        static_cast<cl_float>(velocity.getX()),
        static_cast<cl_float>(velocity.getY())
    };

    m_dataLocal.dirty = true;
}

/* ************************************************************************ */

void Lattice::setOutletDynamics(CoordinateType coord, DensityType density) noexcept
{
    CECE_ASSERT(inRange(coord));

    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());

    Data& data = m_dataLocal.buffer[i];
    data.dynamics = static_cast<cl_int>(Dynamics::Outlet);
    data.outlet.position = coordToBcPosition(getSize(), coord);
    data.outlet.density = density;

    m_dataLocal.dirty = true;
}

/* ************************************************************************ */

void Lattice::setObjectDynamics(CoordinateType coord, VelocityType velocity) noexcept
{
    CECE_ASSERT(inRange(coord));

    const auto i = calcOffset(getSize(), coord.getX(), coord.getY());

    Data& data = m_dataLocal.buffer[i];
    data.dynamics = static_cast<cl_int>(Dynamics::None);

    m_dataLocal.dirty = true;
}

/* ************************************************************************ */

void Lattice::initDefault()
{
    const cl_int2 size = {
        static_cast<cl_int>(getSize().getWidth()),
        static_cast<cl_int>(getSize().getHeight())
    };

    cl_event initEvent;
    cl_event syncEvent;

    {
        const cl_float2 u = {0, 0};
        const cl_float rho = Descriptor::DEFAULT_DENSITY;

        // Set kernel arguments
        clSetKernelArg(m_initKernel, 0, sizeof(size), &size);
        clSetKernelArg(m_initKernel, 1, sizeof(u), &u);
        clSetKernelArg(m_initKernel, 2, sizeof(rho), &rho);
        clSetKernelArg(m_initKernel, 3, sizeof(m_df), &m_df);

        // Enqueue command
        const StaticArray<size_t, 3> dim = {
            getSize().getWidth(),
            getSize().getHeight(),
            Descriptor::SIZE
        };

        CL_CHECK(clEnqueueNDRangeKernel(
            m_commandQueue,
            m_initKernel,
            dim.size(),
            nullptr,
            dim.data(),
            nullptr,
            0, nullptr,
            &initEvent
        ));
    }

    // Calculate density and velocity from new dfs
    {
        CL_CHECK(clSetKernelArg(m_syncKernel, 0, sizeof(size), &size));
        CL_CHECK(clSetKernelArg(m_syncKernel, 1, sizeof(m_df), &m_df));
        CL_CHECK(clSetKernelArg(m_syncKernel, 2, sizeof(m_velocity), &m_velocity));
        CL_CHECK(clSetKernelArg(m_syncKernel, 3, sizeof(m_density), &m_density));
        CL_CHECK(clSetKernelArg(m_syncKernel, 4, sizeof(m_data), &m_data));

        const StaticArray<size_t, 2> dim = {
            getSize().getWidth(),
            getSize().getHeight()
        };

        CL_CHECK(clEnqueueNDRangeKernel(
            m_commandQueue,
            m_syncKernel,
            dim.size(),
            nullptr,
            dim.data(),
            nullptr,
            1, &initEvent,
            &syncEvent
        ));
    }

    // Insert barrier so any following command must wait to initialization
    {
        CL_CHECK(clEnqueueBarrierWithWaitList(
            m_commandQueue,
            1, &syncEvent,
            nullptr
        ));
    }
}

/* ************************************************************************ */

void Lattice::update(unsigned int count)
{
    if (m_dfUpload)
    {
        uploadDf();
        m_dfUpload = false;
    }

    // Update GPU data
    if (m_dataLocal.dirty)
    {
        const size_t size = getSize().getWidth() * getSize().getHeight();

        CL_CHECK(clEnqueueWriteBuffer(
            m_commandQueue,
            m_data,
            CL_FALSE,
            0,
            size * sizeof(Data),
            m_dataLocal.buffer.data(),
            0,
            nullptr,
            nullptr
        ));

        m_dataLocal.dirty = false;
    }

    // Wait until uploades are finished
    CL_CHECK(clEnqueueBarrierWithWaitList(
        m_commandQueue,
        0, nullptr,
        nullptr
    ));

    const cl_int2 size = {
        static_cast<cl_int>(getSize().getWidth()),
        static_cast<cl_int>(getSize().getHeight())
    };

    while (count--)
    {
        cl_event collideEvent;
        cl_event streamEvent;
        cl_event bcEvent;
        cl_event syncEvent;

        // Collide
        {
            const cl_float omega = getOmega();

            CL_CHECK(clSetKernelArg(m_collideKernel, 0, sizeof(size), &size));
            CL_CHECK(clSetKernelArg(m_collideKernel, 1, sizeof(m_velocity), &m_velocity));
            CL_CHECK(clSetKernelArg(m_collideKernel, 2, sizeof(m_density), &m_density));
            CL_CHECK(clSetKernelArg(m_collideKernel, 3, sizeof(m_df), &m_df));
            CL_CHECK(clSetKernelArg(m_collideKernel, 4, sizeof(m_data), &m_data));
            CL_CHECK(clSetKernelArg(m_collideKernel, 5, sizeof(omega), &omega));

            const StaticArray<size_t, 2> dim = {
                getSize().getWidth(),
                getSize().getHeight()
            };

            CL_CHECK(clEnqueueNDRangeKernel(
                m_commandQueue,
                m_collideKernel,
                dim.size(),
                nullptr,
                dim.data(),
                nullptr,
                0, nullptr,
                &collideEvent
            ));
        }

        // Stream
        {
            CL_CHECK(clSetKernelArg(m_streamKernel, 0, sizeof(size), &size));
            CL_CHECK(clSetKernelArg(m_streamKernel, 1, sizeof(m_df), &m_df));
            CL_CHECK(clSetKernelArg(m_streamKernel, 2, sizeof(m_df2), &m_df2));

            const StaticArray<size_t, 3> dim = {
                getSize().getWidth(),
                getSize().getHeight(),
                Descriptor::SIZE
            };

            CL_CHECK(clEnqueueNDRangeKernel(
                m_commandQueue,
                m_streamKernel,
                dim.size(),
                nullptr,
                dim.data(),
                nullptr,
                1, &collideEvent,
                &streamEvent
            ));
        }

        // Swap df buffers
        std::swap(m_df, m_df2);

        // BC
        {
            CL_CHECK(clSetKernelArg(m_bcKernel, 0, sizeof(size), &size));
            CL_CHECK(clSetKernelArg(m_bcKernel, 1, sizeof(m_df), &m_df));
            CL_CHECK(clSetKernelArg(m_bcKernel, 2, sizeof(m_data), &m_data));

            const StaticArray<size_t, 2> dim = {
                getSize().getWidth(),
                getSize().getHeight()
            };

            CL_CHECK(clEnqueueNDRangeKernel(
                m_commandQueue,
                m_bcKernel,
                dim.size(),
                nullptr,
                dim.data(),
                nullptr,
                1, &streamEvent,
                &bcEvent
            ));
        }

        // Sync
        {
            CL_CHECK(clSetKernelArg(m_syncKernel, 0, sizeof(size), &size));
            CL_CHECK(clSetKernelArg(m_syncKernel, 1, sizeof(m_df), &m_df));
            CL_CHECK(clSetKernelArg(m_syncKernel, 2, sizeof(m_velocity), &m_velocity));
            CL_CHECK(clSetKernelArg(m_syncKernel, 3, sizeof(m_density), &m_density));
            CL_CHECK(clSetKernelArg(m_syncKernel, 4, sizeof(m_data), &m_data));

            const StaticArray<size_t, 2> dim = {
                getSize().getWidth(),
                getSize().getHeight()
            };

            CL_CHECK(clEnqueueNDRangeKernel(
                m_commandQueue,
                m_syncKernel,
                dim.size(),
                nullptr,
                dim.data(),
                nullptr,
                1, &bcEvent,
                &syncEvent
            ));
        }

        CL_CHECK(clEnqueueBarrierWithWaitList(
            m_commandQueue,
            1, &syncEvent,
            nullptr
        ));
    }

    // Wait to finish
    CL_CHECK(clFinish(m_commandQueue));
}

/* ************************************************************************ */

void Lattice::uploadDf()
{
    const size_t size = Descriptor::SIZE * getSize().getWidth() * getSize().getHeight();

    cl_event uploadEvent;

    // Write local DF data
    {
        CL_CHECK(clEnqueueWriteBuffer(
            m_commandQueue,
            m_df,
            CL_FALSE,
            0,
            size * sizeof(cl_float),
            m_dfLocal.buffer.data(),
            0,
            nullptr,
            &uploadEvent
        ));
    }

    // Calculate density and velocity from new dfs
    {
        CL_CHECK(clSetKernelArg(m_syncKernel, 0, sizeof(size), &size));
        CL_CHECK(clSetKernelArg(m_syncKernel, 1, sizeof(m_df), &m_df));
        CL_CHECK(clSetKernelArg(m_syncKernel, 2, sizeof(m_velocity), &m_velocity));
        CL_CHECK(clSetKernelArg(m_syncKernel, 3, sizeof(m_density), &m_density));
        CL_CHECK(clSetKernelArg(m_syncKernel, 4, sizeof(m_data), &m_data));

        const StaticArray<size_t, 2> dim = {
            getSize().getWidth(),
            getSize().getHeight()
        };

        CL_CHECK(clEnqueueNDRangeKernel(
            m_commandQueue,
            m_syncKernel,
            dim.size(),
            nullptr,
            dim.data(),
            nullptr,
            1, &uploadEvent,
            nullptr
        ));
    }
}

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
