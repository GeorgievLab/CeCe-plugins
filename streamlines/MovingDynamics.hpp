/* ************************************************************************ */
/* Georgiev Lab (c) 2015-2016                                               */
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
#include "cece/core/ViewPtr.hpp"

// CeCe Plugins
#include "Dynamics.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

/**
 * @brief Moving obstacle dynamics.
 */
class MovingDynamics : public Dynamics
{

// Public Ctors & Dtors
public:


    /**
     * @brief      Constructor.
     *
     * @param      velocity  The obstacle velocity.
     */
    explicit MovingDynamics(VelocityType velocity)
        : m_velocity(velocity)
    {
        // Nothing to do
    }


// Public Operations
public:


    /**
     * @brief Compute macroscopic density in node.
     *
     * @param data Lattice data.
     *
     * @return Macroscopic density.
     */
    DensityType computeDensity(const DataType& data) const noexcept;


    /**
     * @brief Compute macroscopic velocity in node.
     *
     * @param data Lattice data.
     *
     * @return Macroscopic velocity.
     */
    VelocityType computeVelocity(const DataType& data) const noexcept;


    /**
     * @brief Compute macroscopic momentum in node.
     *
     * @param data Lattice data.
     *
     * @return Macroscopic momentum.
     */
    MomentumType computeMomentum(const DataType& data) const noexcept;


    /**
     * @brief Initialize node equilibrum.
     *
     * @param density  Macroscopic density.
     * @param iPop     Direction index.
     * @param velocity Macroscopic velocity.
     *
     * @return Equilibrum distribution.
     */
    DensityType computeEquilibrium(DirectionType iPop, DensityType density,
        VelocityType velocity) const noexcept;


    /**
     * @brief Initialize node equilibrum.
     *
     * @param data     Lattice data.
     * @param velocity Macroscopic velocity.
     * @param density  Macroscopic density.
     */
    void initEquilibrium(DataType& data, VelocityType velocity, DensityType density) const noexcept;


    /**
     * @brief Define node macroscopic density.
     *
     * @param data    Lattice data.
     * @param density Required macroscopic density.
     */
    void defineDensity(DataType& data, DensityType density) const noexcept;


    /**
     * @brief Define node macroscopic velocity.
     *
     * @param data    Lattice data.
     * @param velocity Required macroscopic velocity.
     */
    void defineVelocity(DataType& data, VelocityType velocity) const noexcept;


    /**
     * @brief Perform node collision.
     *
     * @param data Lattice data.
     */
    void collide(DataType& data) const noexcept;


    /**
     * @brief      Prepare pool of dynamics.
     *
     * @param[in]  size  The maximum size.
     */
    static void poolPrepare(int size) noexcept;


    /**
     * @brief      Create a instance of moving dynamics. It's allocated from
     *             pool.
     *
     * @param[in]  velocity  The velocity.
     *
     * @return     Pointer to dynamics.
     */
    static ViewPtr<MovingDynamics> poolCreate(VelocityType velocity) noexcept;


    /**
     * @brief      Clear pool of dynamics.
     */
    static void poolClear() noexcept;


// Private Data Members
private:

    /// Obstacle velocity.
    VelocityType m_velocity;

};

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
