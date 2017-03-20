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

// C++
#include <cstddef>
#include <utility>

// CeCe
#include "cece/core/VectorRange.hpp"

// Plugin
#include "../Descriptor.hpp"

/* ************************************************************************ */

namespace cece {
namespace plugin {
namespace streamlines {
namespace cpu {

/* ************************************************************************ */

Lattice::Lattice(SizeType size, RealType omega)
    : streamlines::Lattice(size, omega)
    , m_data(size)
{
    // Nothing to do
}

/* ************************************************************************ */

Dynamics Lattice::getDynamics(CoordinateType coord) const noexcept
{
    return m_data[coord].getDynamics();
}

/* ************************************************************************ */

Lattice::VelocityType Lattice::getVelocity(CoordinateType coord) const noexcept
{
    return m_data[coord].computeVelocity();
}

/* ************************************************************************ */

void Lattice::setVelocity(CoordinateType coord, VelocityType velocity) noexcept
{
    return m_data[coord].defineVelocity(velocity);
}

/* ************************************************************************ */

Lattice::DensityType Lattice::getDensity(CoordinateType coord) const noexcept
{
    return m_data[coord].computeDensity();
}

/* ************************************************************************ */

void Lattice::setDensity(CoordinateType coord, DensityType density) noexcept
{
    return m_data[coord].defineDensity(density);
}

/* ************************************************************************ */

void Lattice::setVelocityDensity(CoordinateType coord, VelocityType velocity, DensityType density) noexcept
{
    return m_data[coord].initEquilibrium(velocity, density);
}

/* ************************************************************************ */

Lattice::DistributionsType Lattice::getDistributions(CoordinateType coord) const noexcept
{
    return m_data[coord].getDistributions();
}

/* ************************************************************************ */

void Lattice::setDistributions(CoordinateType coord, DistributionsType distributions) noexcept
{
    m_data[coord].setDistributions(distributions);
}

/* ************************************************************************ */

void Lattice::setNoneDynamics(CoordinateType coord) noexcept
{
    m_data[coord].setNoneDynamics();
}

/* ************************************************************************ */

void Lattice::setFluidDynamics(CoordinateType coord) noexcept
{
    m_data[coord].setFluidDynamics();
}

/* ************************************************************************ */

void Lattice::setWallDynamics(CoordinateType coord) noexcept
{
    m_data[coord].setWallDynamics();
}

/* ************************************************************************ */

void Lattice::setInletDynamics(CoordinateType coord, VelocityType velocity) noexcept
{
    if (coord.getX() == 0)
        m_data[coord].setInletDynamics(Node::BoundaryPosition::Left, velocity);
    else if (coord.getX() == getSize().getX() - 1)
        m_data[coord].setInletDynamics(Node::BoundaryPosition::Right, velocity);
    else if (coord.getY() == 0)
        m_data[coord].setInletDynamics(Node::BoundaryPosition::Top, velocity);
    else if (coord.getY() == getSize().getY() - 1)
        m_data[coord].setInletDynamics(Node::BoundaryPosition::Bottom, velocity);
}

/* ************************************************************************ */

void Lattice::setOutletDynamics(CoordinateType coord, DensityType density) noexcept
{
    if (coord.getX() == 0)
        m_data[coord].setOutletDynamics(Node::BoundaryPosition::Left, density);
    else if (coord.getX() == getSize().getX() - 1)
        m_data[coord].setOutletDynamics(Node::BoundaryPosition::Right, density);
    else if (coord.getY() == 0)
        m_data[coord].setOutletDynamics(Node::BoundaryPosition::Top, density);
    else if (coord.getY() == getSize().getY() - 1)
        m_data[coord].setOutletDynamics(Node::BoundaryPosition::Bottom, density);
}

/* ************************************************************************ */

void Lattice::setObjectDynamics(CoordinateType coord, VelocityType velocity) noexcept
{
    m_data[coord].setObjectDynamics(velocity);
}

/* ************************************************************************ */

void Lattice::initDefault()
{
    for (auto& node : m_data)
        node.initEquilibrium();
}

/* ************************************************************************ */

void Lattice::update(unsigned int count)
{
    removeUnreachableDynamics(Dynamics::Wall);

    while (count--)
        collideAndStream();
}

/* ************************************************************************ */

void Lattice::collide()
{
    constexpr Descriptor::PopIndexType half = (Descriptor::SIZE - 1) / 2;

    for (auto&& c : range(getSize()))
    {
        auto& cell = get(c);
        cell.collide(getOmega());

        // TODO: block swap
        for (Descriptor::PopIndexType i = 1; i <= half; ++i)
        {
            using std::swap;
            swap(cell[i], cell[i + half]);
        }
    }
}

/* ************************************************************************ */

void Lattice::stream()
{
    constexpr Descriptor::PopIndexType half = (Descriptor::SIZE - 1) / 2;

    for (auto&& c : range(getSize()))
    {
        for (Descriptor::PopIndexType i = 1; i <= half; ++i)
        {
            // Calculate new coordinates
            const auto newCoord = c + Vector<unsigned int>(Descriptor::DIRECTION_VELOCITIES[i]);

            // Swap
            if (inRange(newCoord))
            {
                using std::swap;
                swap(get(c)[i + half], get(newCoord)[i]);
            }
        }
    }
}

/* ************************************************************************ */

void Lattice::collideAndStream()
{
    collide();
    stream();
/*
    // FIXME: something's wrong
    constexpr Descriptor::IndexType half = (Descriptor::SIZE - 1) / 2;
    const auto size = getSize();

    //for (auto&& c : range(getSize()))
    for (CoordinateType::ValueType x = 0; x < size.getWidth(); ++x)
    for (CoordinateType::ValueType y = 0; y < size.getHeight(); ++y)
    {
        const CoordinateType c{x, y};

        auto& cell = get(c);
        cell.collide();

        for (Descriptor::IndexType i = 1; i <= half; ++i)
        {
            // Calculate new coordinates
            const Vector<Descriptor::IndexType> newCoord = c + Descriptor::DIRECTION_VELOCITIES[i];

            if (inRange(newCoord))
            {
                auto& cellNext = get(newCoord);

                const auto tmp = cell[i];
                cell[i] = cell[i + half];
                cell[i + half] = cellNext[i];
                cellNext[i] = tmp;
            }
        }
    }
*/
}

/* ************************************************************************ */

void Lattice::removeUnreachableDynamics(Dynamics dynamics)
{
    using Offset = Vector<typename std::make_signed<Descriptor::PopIndexType>::type>;

    static const StaticArray<Offset, 9> OFFSETS{{
        Offset{ 0,  0},
        Offset{ 1,  0}, Offset{-1,  0}, Offset{ 0,  1}, Offset{ 1,  1},
        Offset{-1,  1}, Offset{ 0, -1}, Offset{ 1, -1}, Offset{-1, -1}
    }};

    // Foreach all cells
    for (auto&& c : range(getSize()))
    {
        if (getDynamics(c) != dynamics)
            continue;

        bool test = true;

        for (std::size_t i = 0; i < OFFSETS.size(); ++i)
        {
            // TODO: simplify
            Dynamics type = Dynamics::None;

            if (inRange(c + OFFSETS[i]))
                type = getDynamics(c + OFFSETS[i]);

            test = test && (
                type == Dynamics::None ||
                type == dynamics
            );
        }

        if (test)
            setNoneDynamics(c);
    }
}

/* ************************************************************************ */

}
}
}
}

/* ************************************************************************ */
