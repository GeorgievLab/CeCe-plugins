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

// Declaration
#include "MovingDynamics.hpp"

// C++
#include <utility>
#include <iterator>
#include <numeric>

// CeCe
#include "cece/core/Zero.hpp"

// Plugin
#include "config.hpp"
#include "Descriptor.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

DynamicArray<MovingDynamics> g_pool;

/* ************************************************************************ */

}

/* ************************************************************************ */

MovingDynamics::DensityType
MovingDynamics::computeDensity(const DataType& data) const noexcept
{
    using std::begin;
    using std::end;
    return std::accumulate(begin(data), end(data), DensityType(0.0));
}

/* ************************************************************************ */

MovingDynamics::VelocityType
MovingDynamics::computeVelocity(const DataType& data) const noexcept
{
    return Zero;
}

/* ************************************************************************ */

MovingDynamics::MomentumType
MovingDynamics::computeMomentum(const DataType& data) const noexcept
{
    return Zero;
}

/* ************************************************************************ */

MovingDynamics::DensityType
MovingDynamics::computeEquilibrium(DirectionType iPop, DensityType density,
    VelocityType velocity) const noexcept
{
    return Descriptor::calcEquilibrium(
        Descriptor::DIRECTION_WEIGHTS[iPop],
        Descriptor::DIRECTION_VELOCITIES[iPop],
        density,
        velocity
    );
}

/* ************************************************************************ */

void
MovingDynamics::initEquilibrium(DataType& data, VelocityType velocity, DensityType density) const noexcept
{
    for (Descriptor::DirectionType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
        data[iPop] = computeEquilibrium(iPop, density, velocity);
}

/* ************************************************************************ */

void
MovingDynamics::defineDensity(DataType& data, DensityType density) const noexcept
{
    // Nothing to do
}

/* ************************************************************************ */

void
MovingDynamics::defineVelocity(DataType& data, VelocityType velocity) const noexcept
{
    // Nothing to do
}

/* ************************************************************************ */

void
MovingDynamics::collide(DataType& data) const noexcept
{
    DataType temp;

    // Move updated values into opposite directions
    for (Descriptor::DirectionType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        const auto wi = Descriptor::DIRECTION_WEIGHTS[iPop];
        const auto ui = Descriptor::DIRECTION_VELOCITIES[iPop];
        const auto iop = Descriptor::DIRECTION_OPPOSITES[iPop];

        temp[iPop] = data[iop] + 6 * wi * dot(ui, m_velocity);

        // Result value must be positive
        CECE_ASSERT(temp[iPop] > 0);
    }

    // Copy updated values
    data = temp;
}

/* ************************************************************************ */

void MovingDynamics::poolPrepare(int size) noexcept
{
    g_pool.reserve(size);
}

/* ************************************************************************ */

ViewPtr<MovingDynamics> MovingDynamics::poolCreate(VelocityType velocity) noexcept
{
    g_pool.emplace_back(velocity);

    return &g_pool.back();
}

/* ************************************************************************ */

void MovingDynamics::poolClear() noexcept
{
    g_pool.clear();
}

/* ************************************************************************ */

}
}
}

/* ************************************************************************ */
