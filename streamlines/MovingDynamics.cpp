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

const StaticArray<StaticArray<Descriptor::DirectionType, 3>, 3> VEC_TO_INDEX{{
    {{3, 4, 5}},
    {{2, 0, 6}},
    {{1, 8, 7}}
}};

/* ************************************************************************ */

const StaticArray<StaticArray<Descriptor::DirectionType, Descriptor::SIZE>, Descriptor::SIZE> NORMAL_MAPPING{{
    {{0, 1, 2, 3, 4, 5, 6, 7, 8}},
    {{0, 5, 4, 3, 4, 5, 6, 7, 6}},
    {{0, 7, 6, 5, 4, 5, 6, 7, 8}},
    {{0, 1, 8, 7, 6, 5, 6, 7, 8}},
    {{0, 1, 2, 1, 8, 7, 6, 7, 8}},
    {{0, 1, 2, 3, 2, 1, 8, 7, 8}},
    {{0, 1, 2, 3, 4, 3, 2, 1, 8}},
    {{0, 1, 2, 3, 4, 5, 4, 3, 2}},
    {{0, 3, 2, 3, 4, 5, 6, 5, 4}},
}};

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
    const auto n = Vector<int>(
        std::lround(m_normal.getX()),
        std::lround(m_normal.getY())
    );

    // Normal index
    const auto ni = VEC_TO_INDEX[n.getY() + 1][n.getX() + 1];

    DataType temp;

    RealType wp = 1.0;
    RealType wr = 0;

    // Move updated values into opposite directions
    for (Descriptor::DirectionType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        const auto wi = Descriptor::DIRECTION_WEIGHTS[iPop];
        const auto ui = Descriptor::DIRECTION_VELOCITIES[iPop];
        const auto iop = Descriptor::DIRECTION_OPPOSITES[iPop];
        const auto ir = NORMAL_MAPPING[ni][iPop];

        temp[iPop] =
            wp * (data[iop] + 2 * wi * 3 * dot(ui, m_velocity)) +
            (1.0 - wp) * (data[ir] + wr * 2 * wi * 3 * dot(ui, m_normal * dot(m_normal, m_normal)))
        ;

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

ViewPtr<MovingDynamics> MovingDynamics::poolCreate(VelocityType velocity, Vector<RealType> normal) noexcept
{
    g_pool.emplace_back(velocity, normal);

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
