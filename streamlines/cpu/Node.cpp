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
#include "Node.hpp"

// C++
#include <iterator>
#include <numeric>

// CeCe
#include "cece/core/Map.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace cpu {

/* ************************************************************************ */

namespace {

/* ************************************************************************ */

const Map<Node::BoundaryPosition, StaticArray<Descriptor::PopIndexType, 3>> CENTER_RHO{{
    {Node::BoundaryPosition::Right,    Descriptor::MIDDLE_COLUMN},
    {Node::BoundaryPosition::Left,     Descriptor::MIDDLE_COLUMN},
    {Node::BoundaryPosition::Top,      Descriptor::MIDDLE_LINE},
    {Node::BoundaryPosition::Bottom,   Descriptor::MIDDLE_LINE}
}};

/* ************************************************************************ */

const Map<Node::BoundaryPosition, StaticArray<Descriptor::PopIndexType, 3>> KNOWN_RHO{{
    {Node::BoundaryPosition::Right,    Descriptor::RIGHT_COLUMN},
    {Node::BoundaryPosition::Left,     Descriptor::LEFT_COLUMN},
    {Node::BoundaryPosition::Top,      Descriptor::TOP_LINE},
    {Node::BoundaryPosition::Bottom,   Descriptor::BOTTOM_LINE}
}};

/* ************************************************************************ */

const Map<Node::BoundaryPosition, StaticArray<Descriptor::PopIndexType, 3>> UNKNOWN_RHO{{
    {Node::BoundaryPosition::Right,    Descriptor::LEFT_COLUMN},
    {Node::BoundaryPosition::Left,     Descriptor::RIGHT_COLUMN},
    {Node::BoundaryPosition::Top,      Descriptor::BOTTOM_LINE},
    {Node::BoundaryPosition::Bottom,   Descriptor::TOP_LINE}
}};

/* ************************************************************************ */

const Map<Node::BoundaryPosition, Node::VelocityType> VELOCITIES{{
    {Node::BoundaryPosition::Right,    {-1,  0}},
    {Node::BoundaryPosition::Left,     { 1,  0}},
    {Node::BoundaryPosition::Top,      { 0, -1}},
    {Node::BoundaryPosition::Bottom,   { 0,  1}}
}};

/* ************************************************************************ */

const Map<Node::BoundaryPosition, Descriptor::PopIndexType> BC_CENTER{{
    {Node::BoundaryPosition::Right,    Descriptor::INDEX_MAP[1][0]},
    {Node::BoundaryPosition::Left,     Descriptor::INDEX_MAP[1][2]},
    {Node::BoundaryPosition::Top,      Descriptor::INDEX_MAP[2][1]},
    {Node::BoundaryPosition::Bottom,   Descriptor::INDEX_MAP[0][1]}
}};

/* ************************************************************************ */

const Map<Node::BoundaryPosition, StaticArray<Descriptor::PopIndexType, 2>> BC_SIDE1{{
    {Node::BoundaryPosition::Right,    {{Descriptor::INDEX_MAP[0][0], Descriptor::INDEX_MAP[2][1]}}},
    {Node::BoundaryPosition::Left,     {{Descriptor::INDEX_MAP[0][2], Descriptor::INDEX_MAP[2][1]}}},
    {Node::BoundaryPosition::Top,      {{Descriptor::INDEX_MAP[2][2], Descriptor::INDEX_MAP[1][0]}}},
    {Node::BoundaryPosition::Bottom,   {{Descriptor::INDEX_MAP[0][2], Descriptor::INDEX_MAP[1][0]}}}
}};

/* ************************************************************************ */

const Map<Node::BoundaryPosition, StaticArray<Descriptor::PopIndexType, 2>> BC_SIDE2{{
    {Node::BoundaryPosition::Right,    {{Descriptor::INDEX_MAP[2][0], Descriptor::INDEX_MAP[0][1]}}},
    {Node::BoundaryPosition::Left,     {{Descriptor::INDEX_MAP[2][2], Descriptor::INDEX_MAP[0][1]}}},
    {Node::BoundaryPosition::Top,      {{Descriptor::INDEX_MAP[2][0], Descriptor::INDEX_MAP[1][2]}}},
    {Node::BoundaryPosition::Bottom,   {{Descriptor::INDEX_MAP[0][0], Descriptor::INDEX_MAP[1][2]}}}
}};

/* ************************************************************************ */

/**
 * @brief      Compute of total sum of given value indices.
 *
 * @param      df    The data.
 * @param      list  List of indices to sum.
 *
 * @tparam     N     Number of indices.
 *
 * @return     Sum of values.
 */
template<size_t N>
Node::DensityType sumValues(Node::DistributionsType& df, StaticArray<Descriptor::PopIndexType, N> list) noexcept
{
    Node::DensityType density{};

    for (auto iPop : list)
        density += df[iPop];

    return density;
}

/* ************************************************************************ */

void zouHeInit(Node::BoundaryPosition position, Node::DistributionsType& df, Node::VelocityType velocity, Node::DensityType density) noexcept
{
    const auto center = BC_CENTER.at(position);
    const auto side1 = BC_SIDE1.at(position);
    const auto side2 = BC_SIDE2.at(position);

    auto eqDiff = [&density, &velocity] (Descriptor::PopIndexType iPop) {
        return
            Descriptor::calcEquilibrium(
                Descriptor::DIRECTION_WEIGHTS[iPop],
                Descriptor::DIRECTION_VELOCITIES[iPop],
                density,
                velocity
            ) - Descriptor::calcEquilibrium(
                Descriptor::DIRECTION_WEIGHTS[Descriptor::opposite(iPop)],
                Descriptor::DIRECTION_VELOCITIES[Descriptor::opposite(iPop)],
                density,
                velocity
            )
        ;
    };

    auto dataDiff = [&df] (Descriptor::PopIndexType iPop) {
        return df[iPop] - df[Descriptor::opposite(iPop)];
    };

    const auto side1_0 = side1[0];
    const auto side1_1 = side1[1];
    const auto side2_0 = side2[0];
    const auto side2_1 = side2[1];

    // Center
    df[center] = df[Descriptor::opposite(center)]
        + eqDiff(center)
    ;
    CECE_ASSERT(df[center] > 0);

    // Side 1
    df[side1_0] = df[Descriptor::opposite(side1_0)]
        + eqDiff(side1_0)
        + 0.5 * dataDiff(side1_1)
    ;
    CECE_ASSERT(df[side1_0] > 0);

    // Side 2
    df[side2_0] = df[Descriptor::opposite(side2_0)]
        + eqDiff(side2_0)
        + 0.5 * dataDiff(side2_1)
    ;
    CECE_ASSERT(df[side2_0] > 0);
}

/* ************************************************************************ */

}

/* ************************************************************************ */

void Node::setNoneDynamics() noexcept
{
    m_dynamics = Dynamics::None;
}

/* ************************************************************************ */

void Node::setFluidDynamics() noexcept
{
    m_dynamics = Dynamics::Fluid;
}

/* ************************************************************************ */

void Node::setWallDynamics() noexcept
{
    m_dynamics = Dynamics::Wall;
}

/* ************************************************************************ */

void Node::setInletDynamics(BoundaryPosition position, VelocityType velocity) noexcept
{
    m_dynamics = Dynamics::Inlet;
    m_data.inlet.position = position;
    m_data.inlet.velocity = velocity;
}

/* ************************************************************************ */

void Node::setOutletDynamics(BoundaryPosition position, DensityType density) noexcept
{
    m_dynamics = Dynamics::Outlet;
    m_data.outlet.position = position;
    m_data.outlet.density = density;
}

/* ************************************************************************ */

void Node::setObjectDynamics(VelocityType velocity) noexcept
{
    m_dynamics = Dynamics::Object;
    m_data.object.velocity = velocity;
}

/* ************************************************************************ */

void Node::defineVelocity(const VelocityType& velocity) noexcept
{
    initEquilibrium(velocity, computeDensity());
}

/* ************************************************************************ */

void Node::defineDensity(DensityType density) noexcept
{
    initEquilibrium(computeVelocity(), density);
}

/* ************************************************************************ */

Node::DensityType Node::computeDensity() const noexcept
{
    using std::begin;
    using std::end;
    return std::accumulate(begin(m_df), end(m_df), DensityType(0.0));
}

/* ************************************************************************ */

Node::MomentumType Node::computeMomentum() const noexcept
{
    MomentumType momentum = Zero;

    for (PopIndexType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
        momentum += Descriptor::DIRECTION_VELOCITIES[iPop] * m_df[iPop];

    return momentum;
}

/* ************************************************************************ */

Node::VelocityType Node::computeVelocity() const noexcept
{
    return computeMomentum() / computeDensity();
}

/* ************************************************************************ */

void Node::initEquilibrium(VelocityType velocity, DensityType density) noexcept
{
    for (PopIndexType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        m_df[iPop] = Descriptor::calcEquilibrium(
            Descriptor::DIRECTION_WEIGHTS[iPop],
            Descriptor::DIRECTION_VELOCITIES[iPop],
            density,
            velocity
        );

        CECE_ASSERT(m_df[iPop] > 0);
    }
}

/* ************************************************************************ */

void Node::collide(RealType omega) noexcept
{
    switch (getDynamics())
    {
    case Dynamics::Inlet:
        boundaryInlet();
        break;

    case Dynamics::Outlet:
        boundaryOutlet();
        break;

    default:
        break;
    }

    switch (getDynamics())
    {
    case Dynamics::Fluid:
    case Dynamics::Inlet:
    case Dynamics::Outlet:
        collideFluid(omega);
        break;

    case Dynamics::Object:
        collideObject(omega);
        break;

    case Dynamics::Wall:
        collideWall();
        break;

    case Dynamics::None:
        // Nothing to do
        break;
    }
}

/* ************************************************************************ */

void Node::collideFluid(RealType omega) noexcept
{
    const auto density = computeDensity();
    CECE_ASSERT(density > 0);
    const auto velocity = computeVelocity();

    // Perform collision for all populations
    for (PopIndexType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        // Calculate equilibrium distribution
        const auto feq = Descriptor::calcEquilibrium(
            Descriptor::DIRECTION_WEIGHTS[iPop],
            Descriptor::DIRECTION_VELOCITIES[iPop],
            density,
            velocity
        );

        CECE_ASSERT(feq > 0);
        CECE_ASSERT(m_df[iPop] > 0);

        // Collide
        m_df[iPop] += -omega * (m_df[iPop] - feq);

        // Result value must be positive
        CECE_ASSERT(m_df[iPop] > 0);
    }
}

/* ************************************************************************ */

void Node::collideWall() noexcept
{
    constexpr PopIndexType half = (Descriptor::SIZE - 1) / 2;

    // Swap opposite (without temporary memory)
    for (PopIndexType i = 1; i <= half; ++i)
    {
        using std::swap;
        swap(m_df[i], m_df[i + half]);
    }
}

/* ************************************************************************ */

void Node::collideObject(RealType omega) noexcept
{
    DistributionsType temp;

    // Move updated values into opposite directions
    for (Descriptor::PopIndexType iPop = 0; iPop < Descriptor::SIZE; ++iPop)
    {
        const auto wi = Descriptor::DIRECTION_WEIGHTS[iPop];
        const auto ui = Descriptor::DIRECTION_VELOCITIES[iPop];
        const auto iop = Descriptor::DIRECTION_OPPOSITES[iPop];

        temp[iPop] = m_df[iop] + 6 * wi * dot(ui, m_data.object.velocity);

        // Result value must be positive
        CECE_ASSERT(temp[iPop] > 0);
    }

    // Copy updated values
    m_df = temp;
}

/* ************************************************************************ */

void Node::boundaryInlet() noexcept
{
    const auto center = sumValues(m_df, CENTER_RHO.at(m_data.inlet.position));
    const auto known  = sumValues(m_df, KNOWN_RHO.at(m_data.inlet.position));
    const auto velP = m_data.inlet.velocity.dot(VELOCITIES.at(m_data.inlet.position));

    const auto density = 1.0 / (1.0 - velP) * (center + 2 * known);

    zouHeInit(m_data.inlet.position, m_df, m_data.inlet.velocity, density);
}

/* ************************************************************************ */

void Node::boundaryOutlet() noexcept
{
    const auto center = sumValues(m_df, CENTER_RHO.at(m_data.outlet.position));
    const auto known  = sumValues(m_df, KNOWN_RHO.at(m_data.outlet.position));

    // Speed
    const RealType speed = (1.0 - 1.0 / m_data.outlet.density * (center + 2 * known));

    // Velocity vector
    const VelocityType velocity = speed * VELOCITIES.at(m_data.outlet.position);

    zouHeInit(m_data.outlet.position, m_df, velocity, m_data.outlet.density);
}

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
