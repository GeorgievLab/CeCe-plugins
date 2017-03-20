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

// CeCe
#include "cece/core/DynamicArray.hpp"

// Plugin
#include "../Lattice.hpp"

// OpenCL
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace gpu {

/* ************************************************************************ */

/**
 * @brief      GPU implementation of Lattice Boltzman.
 */
class Lattice : public streamlines::Lattice
{

// Public Enums
public:


    /**
     * @brief      Boundary position.
     */
    enum class BoundaryPosition
    {
        Right  = 0,
        Left   = 1,
        Top    = 2,
        Bottom = 3
    };

// Public Types
public:


    /// Real type used for computation.
    using cl_RealType = cl_double;

    /// Real type 2D vector used for computation.
    using cl_RealType2 = cl_double2;


// Public Ctors & Dtors
public:


    /**
     * @brief      Constructor.
     *
     * @param[in]  size   The lattice size.
     * @param[in]  omega  The relaxation frequency.
     */
    Lattice(SizeType size, RealType omega);


    /**
     * @brief      Destructor.
     */
    ~Lattice();


// Public Accessors & Mutators
public:


    /**
     * @brief      Get dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The dynamics.
     */
    Dynamics getDynamics(CoordinateType coord) const noexcept override;


    /**
     * @brief      Get velocity at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The velocity.
     */
    VelocityType getVelocity(CoordinateType coord) const noexcept override;


    /**
     * @brief      Set velocity at given coordinate.
     *
     * @details    Function update node velocity using equilibrium distribution
     *             with node's density.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The velocity.
     */
    void setVelocity(CoordinateType coord, VelocityType velocity) noexcept override;


    /**
     * @brief      Get density at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     *
     * @return     The density.
     */
    DensityType getDensity(CoordinateType coord) const noexcept override;


    /**
     * @brief      Set density at given coordinate.
     *
     * @details    Function update node density using equilibrium distribution
     *             with node's velocity.
     *
     * @param[in]  coord    The coordinate.
     * @param[in]  density  The density.
     */
    void setDensity(CoordinateType coord, DensityType density) noexcept override;


    /**
     * @brief      Set velocity and density at given coordinate.
     *
     * @details    Function update node velocity and density using equilibrium
     *             distribution.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The velocity.
     * @param[in]  density   The density.
     */
    void setVelocityDensity(CoordinateType coord, VelocityType velocity, DensityType density) noexcept override;


    /**
     * @brief      Get distribution functions at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     *
     * @return     The distribution functions.
     */
    DistributionsType getDistributions(CoordinateType coord) const noexcept override;


    /**
     * @brief      Set distribution functions at given coordinate.
     *
     * @param[in]  coord          The coordinate.
     * @param[in]  distributions  The distribution functions.
     */
    void setDistributions(CoordinateType coord, DistributionsType distributions) noexcept override;


    /**
     * @brief      Set none dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    void setNoneDynamics(CoordinateType coord) noexcept override;


    /**
     * @brief      Set fluid dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    void setFluidDynamics(CoordinateType coord) noexcept override;


    /**
     * @brief      Set wall dynamics at given coordinate.
     *
     * @param[in]  coord  The coordinate.
     */
    void setWallDynamics(CoordinateType coord) noexcept override;


    /**
     * @brief      Set inlet dynamics at given coordinate.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The inlet velocity.
     */
    void setInletDynamics(CoordinateType coord, VelocityType velocity) noexcept override;


    /**
     * @brief      Set outlet dynamics at given coordinate.
     *
     * @param[in]  coord    The coordinate.
     * @param[in]  density  The outlet density.
     */
    void setOutletDynamics(CoordinateType coord, DensityType density) noexcept override;


    /**
     * @brief      Set an dynamics at given coordinate object.
     *
     * @param[in]  coord     The coordinate.
     * @param[in]  velocity  The object velocity.
     */
    void setObjectDynamics(CoordinateType coord, VelocityType velocity) noexcept override;


// Public Operations
public:


    /**
     * @brief      Initialize whole lattice with default values.
     *
     * @details    As default values it means zero velocity with default density
     *             (equilibrium).
     */
    void initDefault() override;


    /**
     * @brief      Update lattice according to set parameters and current state.
     *
     * @param[in]  count  The number of inner iterations.
     */
    void update(unsigned int count = 1) override;


// Private Structures
private:


    /**
     * @brief      Local buffer structure.
     *
     * @tparam     T     Element type.
     */
    template<typename T>
    struct LocalBuffer
    {
        DynamicArray<T> buffer;
        bool dirty = true;
    };


    /// Fluid dynamics data
    struct FluidData { };

    /// Wall dynamics data
    struct WallData { };

    /// Inlet dynamics data
    struct InletData
    {
        cl_int position;
        cl_RealType2 velocity;
    };

    /// Outlet dynamics data
    struct OutletData
    {
        cl_int position;
        cl_RealType density;
    };

    /// Object dynamics data
    struct ObjectData
    {
        cl_RealType2 velocity;
    };

    /// Dynamics data
    struct Data
    {
        cl_int dynamics;
        struct FluidData fluid;
        struct WallData wall;
        struct InletData inlet;
        struct OutletData outlet;
        struct ObjectData object;
    };

// Private Operations
private:


    /**
     * @brief      Upload distribution functions.
     */
    void uploadDf();


// Private Data Members
private:

    /// OpenCL platform.
    cl_platform_id m_platform = nullptr;

    /// OpenCL device.
    cl_device_id m_device = nullptr;

    /// OpenCL context.
    cl_context m_context = nullptr;

    /// Command queue.
    cl_command_queue m_commandQueue = nullptr;

    /// OpenCL LB program.
    cl_program m_program = nullptr;

    /// Init kernel.
    cl_kernel m_initKernel = nullptr;

    /// Collide kernel.
    cl_kernel m_collideKernel = nullptr;

    /// Stream kernel.
    cl_kernel m_streamKernel = nullptr;

    /// Boundary conditions kernel.
    cl_kernel m_bcKernel = nullptr;

    /// Sync kernel.
    cl_kernel m_syncKernel = nullptr;

    /// Memory for distribution functions.
    cl_mem m_df = nullptr;

    /// Memory for distribution functions.
    cl_mem m_df2 = nullptr;

    /// Memory for velocities.
    cl_mem m_velocity = nullptr;

    /// Memory for densities.
    cl_mem m_density = nullptr;

    /// Memory for data.
    cl_mem m_data = nullptr;

    /// Local distribution functions.
    mutable LocalBuffer<cl_RealType> m_dfLocal;

    /// Local velocities.
    mutable LocalBuffer<cl_RealType2> m_velocityLocal;

    /// Local densities.
    mutable LocalBuffer<cl_RealType> m_densityLocal;

    /// Local data.
    mutable LocalBuffer<Data> m_dataLocal;

    /// If distribution functions should be uploaded to GPU
    mutable bool m_dfUpload = false;
};

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace gpu {

/* ************************************************************************ */

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
